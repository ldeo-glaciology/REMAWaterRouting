import numpy as np
import sklearn.linear_model
import skimage.morphology
import skimage.segmentation
import scipy.ndimage
import dask
import geopandas as gpd
import shapely
from rasterio import RasterioIOError
from tqdm.autonotebook import tqdm
import xarray as xr
import os

#read in the shapefiles of ice shelf grounding lines
IS = gpd.read_file('IceShelf_Antarctica_v02/IceShelf_Antarctica_v02.shp') 
#read in the REMA tile index
REMA_index = gpd.read_file('REMA_Tile_Index/REMA_Tile_Index_Rel1_1.shp')

def download_REMA(shelf):
    #bounding box of ice shelf
    SelectedIceShelf = IS[IS.NAME==shelf]
    [minx,miny,maxx,maxy]= SelectedIceShelf.bounds.values.tolist()[0]

    bbox = shapely.geometry.asPolygon([[minx,miny],[maxx,miny],[maxx,maxy],[minx,maxy],[minx,miny]])

    IS_intersection = np.argwhere(REMA_index.overlaps(bbox).tolist())

    IS_tiles = REMA_index.tile[IS_intersection.flatten()]

    row=np.zeros((len(IS_tiles),1))
    col=np.zeros((len(IS_tiles),1))
    for i in np.arange(0,len(IS_tiles)):
        [row[i],col[i]] = str.split(IS_tiles.to_list()[i],sep='_')

    row = np.int_(row)
    col = np.int_(col)

    #Load the REMA tiles lazily
    uri_fmt = 'https://storage.googleapis.com/pangeo-pgc/8m/{i_idx:02d}_{j_idx:02d}/{i_idx:02d}_{j_idx:02d}_8m_dem_COG_LZW.tif'

    chunksize = 8 * 512
    rows = []
    for i in tqdm(range(row.max()-1, row.min()-1, -1)): #this tile range is Amery Ice Shelf
        cols = []
        for j in range(col.min(),col.max()):
            uri = uri_fmt.format(i_idx=i, j_idx=j)
            try:
                dset = xr.open_rasterio(uri, chunks=chunksize)
                dset_masked = dset.where(dset > 0.0)
                cols.append(dset_masked)
                #print(uri)
            except RasterioIOError:
                pass
        rows.append(cols)

    dsets_rows = [xr.concat(row, 'x') for row in rows]
    ds = xr.concat(dsets_rows, 'y', )
    ds_array = ds.squeeze()
    coarsed_masked_array = ds_array.coarsen(x=100,y=100).mean()
    cellsize = ds_array.res[0]*100*ds_array.res[1]*100
    sample = coarsed_masked_array.chunk((250,250))
    return(sample)

def flag_nans(dem):
    log_nans = np.isnan(dem);
         # handle NaNs

    if np.any(log_nans):
        flag_nans = 1;
    else:
        flag_nans = 0;
    return flag_nans
    
def identifyflats(dem):
    flats=np.zeros(dem.size)

    if flag_nans==True:
        dem = dask.array.nan_to_num(dem, nan=-np.inf) 
    else:
        dem = dem


    nhood = np.ones((3,3))
    # identify flats
    # flats: logical matrix with true where cells don't have lower neighbors

    if len(dem)>1:
        flats = skimage.morphology.erosion(image=dem,selem=nhood) == dem;

        #remove flats at the border
        flats[0:-1,[1, -1]]  = 0;
        flats[[0 -1],0:-1]  = 0;

        flats = skimage.segmentation.clear_border(flats, buffer_size=2)

       # if flag_nans==1:
        #    # remove flat pixels bordering to nans
         #   flats[skimage.morphology.dilation(image=log_nans,selem = nhood)] = 0
    return flats


def identifysills(dem,flats):
    # identify sills
    sills=np.zeros(dem.size)
    nhood = np.ones((3,3))

    log_nans = np.isnan(dem);

    # find sills and set marker
    if len(dem)>1:
        Imr = -np.inf * np.ones(flats.shape);
        Imr[flats.astype(int)] = 0;
        maskeddem = np.multiply(dem,flats)
        Imr = Imr + maskeddem;
        Imr = (skimage.morphology.dilation(image=Imr,selem=nhood) == dem) & np.logical_not(flats);

        if flag_nans==True:
               Imr[log_nans] = 0;

        sills = Imr
    return sills
    
def identifyinteriorbasins(dem):
    log_nans = np.isnan(dem);

    # identify interior basins
    interiorbasins=np.zeros(dem.size)

    if len(dem)>1:
        interiorbasins = skimage.morphology.local_minima(dem);

    if flag_nans==1:
        interiorbasins = np.bitwise_or(interiorbasins,log_nans);
        if interiorbasins.shape[0]>0:
            interiorbasins = skimage.segmentation.clear_border(interiorbasins, buffer_size=0);
        interiorbasins[log_nans] = 0;
    else:
        if interiorbasins.shape[0]>0:
            interiorbasins = skimage.segmentation.clear_border(interiorbasins);

    return interiorbasins
    
import heapq
from imageio import imread


DIR_STRINGS = ["left", "down", "right", "up"]
DIRS = ((-1, 0), (0, -1), (1, 0), (0, 1))

def map_image_to_costs(D, PreSillPixel):
    """
    Read image data and convert it to a marginal cost function,
    a 2D array containing costs for moving through each pixel.
    This cost field forms the input for the weighted distance transform
    zero costs denote exits, infinite costs denote fully impenetrable obstacles.
    In this example, we follow Mercurial standards: obstacles are in black, exits in green,
    accessible space is in white, less accessible space has less white.
    Adapt to your own needs.
    :param image: String of image file or open file descriptor of image
    :return: 2D array representing the cost field
    """
    
    G = (D-np.min(D))/(np.max(D)-np.min(D))**-1
    if len(PreSillPixel)>0:
        for pair in PreSillPixel:
            PSPx = int(pair[0])
            PSPy = int(pair[1])

        #R = abs(D)*255/np.max(abs(D))

        G[PSPx,PSPy] = 1
        #B = abs(D)*255/np.max(abs(D)
    
    data = G*256
    # Exits are present in all green enough places ("G >> R and G")
    exits = np.where(data >= 255 )
    # Obstacles are in black (so at least G and B must be zero)
    obstacles = np.where(data < 1)
    # Convert image to greyscale
    grey_scales = data
    # Boolean index array for places without exits and obstacles
    space = np.ones(grey_scales.shape, dtype=np.bool)
    space[obstacles] = False
    space[exits] = False
    # Cost field: Inversely proportional to greyscale values
    cost_field = np.empty(data.shape)
    cost_field[obstacles] = np.inf
    cost_field[exits] = 0
    cost_field[space] = 1. / (grey_scales[space])
    return cost_field

def _wdt_python(cost_field):
    """
    See `get_weighted_distance_transform`
    :param cost_field: 2D array
    :return: Weighted distance transform array with same shape as `cost_field`
    """
    nx, ny = cost_field.shape
    # Cost for moving along horizontal lines
    costs_x = np.ones([nx + 1, ny], order='F') * np.inf
    costs_x[1:-1, :] = (cost_field[1:, :] + cost_field[:-1, :]) / 2
    # Cost for moving along vertical lines
    costs_y = np.ones([nx, ny + 1], order='F') * np.inf
    costs_y[:, 1:-1] = (cost_field[:, 1:] + cost_field[:, :-1]) / 2

    # Initialize locations (known/unknown/exit/obstacle)
    weighted_distance_transform = np.ones_like(cost_field, order='F') * np.inf
    exit_locs = np.where(cost_field == 0)
    obstacle_locs = np.where(cost_field == np.inf)
    weighted_distance_transform[exit_locs] = 0

    # Initialize Cell structures
    all_cells = {(i, j) for i in range(nx) for j in range(ny)}
    known_cells = {cell for cell in zip(exit_locs[0], exit_locs[1])}
    unknown_cells = all_cells - known_cells - {cell for cell in zip(obstacle_locs[0], obstacle_locs[1])}
    new_candidate_cells = set()
    for cell in known_cells:
        new_candidate_cells |= _get_new_candidate_cells(cell, unknown_cells)
    cand_heap = [(np.inf, cell) for cell in new_candidate_cells]
    # Loop until all unknown cells have a distance value
    if len(cand_heap)>0:
        while True:
            # by repeatedly looping over the new candidate cells
            for cell in new_candidate_cells:
                # Compute a distance for each cell based on its neighbour cells
                distance = _propagate_distance(cell, [costs_x, costs_y], weighted_distance_transform)
                # Store this value in the heap (for fast lookup)
                # Don't check whether we have the distance already in the heap; check on outcome
                heapq.heappush(cand_heap, (distance, cell))
            # See if the heap contains a good value and if so, add it to the field. If not, finish.
            # Since we can store multiple distance values for one cell, we might need to pop a couple of times
            while True:
                min_distance, best_cell = heapq.heappop(cand_heap)
                if weighted_distance_transform[best_cell] == np.inf:
                    # Got a good one: no assigned distance in wdt yet
                    break
                elif min_distance == np.inf:  # No more finite values; done
                    return weighted_distance_transform
            # Good value found, add to the wdt and
            weighted_distance_transform[best_cell] = min_distance
            unknown_cells.remove(best_cell)
            new_candidate_cells = _get_new_candidate_cells(best_cell, unknown_cells)
    else:
        weighted_distance_transform = scipy.ndimage.distance_transform_edt(cost_field)+1
        return weighted_distance_transform
    """
    Checks whether an index exists an array
    :param index: 2D index tuple
    :return: true if lower than tuple, false otherwise
    """
    return (0 <= index[0] < nx) and (0 <= index[1] < ny)
def _exists(index, nx, ny):
    """
    Checks whether an index exists an array
    :param index: 2D index tuple
    :return: true if lower than tuple, false otherwise
    """
    return (0 <= index[0] < nx) and (0 <= index[1] < ny)


def _get_new_candidate_cells(cell, unknown_cells):
    """
    Compute the new candidate cells (cells for which we have no definite distance value yet
    For more information on the algorithm: check fast marching method
    :param cell: tuple of index; a new cell that has been added to the distance field
    :param unknown_cells: set of tuples; all cells still unknown
    :return: Set of new candidate cells for which to compute the distance
    """
    new_candidate_cells = set()
    for direction in DIRS:
        nb_cell = (cell[0] + direction[0], cell[1] + direction[1])
        if nb_cell in unknown_cells:
            new_candidate_cells.add(nb_cell)
    return new_candidate_cells


def _propagate_distance(cell, costs, wdt_field):
    """
    Compute the weighted distance in a cell using costs and distances in other cells
    :param cell: tuple, index of a candidate cell
    :param costs: list of cost arrays in X and Y direction
    :param wdt_field: the weighted distance transform field up until now
    :return: a approximate distance based on the neighbour cells
    """
    nx, ny = wdt_field.shape
    # Find the minimal directions along a grid cell.
    # Assume left and below are best, then overwrite with right and up if they are better
    adjacent_distances = np.ones(4) * np.inf
    pots_from_axis = [0, 0]  # [x direction, y direction]
    costs_from_axis = [np.inf, np.inf]  #
    for i, dir_s in enumerate(DIR_STRINGS):
        # Direction for which we check the cost
        normal = DIRS[i]
        nb_cell = (cell[0] + normal[0], cell[1] + normal[1])
        if not _exists(nb_cell, nx, ny):
            continue
        pot = wdt_field[nb_cell]
        # distance in that neighbour field
        if dir_s == 'left':
            face_index = (nb_cell[0] + 1, nb_cell[1])
        elif dir_s == 'down':
            face_index = (nb_cell[0], nb_cell[1] + 1)
        else:
            face_index = nb_cell
        # Left/right is x, up/down is y
        cost = costs[i % 2][face_index]
        # Proposed cost along this direction
        adjacent_distances[i] = pot + cost
        # If it is cheaper to go from the opposite direction
        if adjacent_distances[i] < adjacent_distances[(i + 2) % 4]:
            pots_from_axis[i % 2] = pot
            costs_from_axis[i % 2] = cost
        hor_pot, ver_pot = pots_from_axis
        hor_cost, ver_cost = costs_from_axis
        # Coefficients of quadratic equation (upwind discretization)
    a = 1. / hor_cost ** 2 + 1. / ver_cost ** 2
    b = -2 * (hor_pot / hor_cost ** 2 + ver_pot / ver_cost ** 2)
    c = (hor_pot / hor_cost) ** 2 + (ver_pot / ver_cost) ** 2 - 1

    D = b ** 2 - 4 * a * c
    # Largest root represents upwind approximation
    x_high = (2 * c) / (-b - math.sqrt(D+0.001)) #edited to prevent divide by zero errors
    return x_high


def drainagebasins(Z,flats,sills,interiorbasins, cellsize):  
    
    #[Iobj,SILLS,IntBasin] = identifyflats(Z);
    Z = Z.data
    Z_ravel = np.ravel(Z)
    nrc = Z_ravel.shape[0]
    
    Iobj  = flats
    SILLS = sills
    IntBasin = interiorbasins

    # Here we choose the distance transform from outside the lakes to the inside and take the locations as sills where the distance is maximum.
    DD = scipy.ndimage.distance_transform_edt(np.bitwise_not(IntBasin));
    MaxIntIX = [0,0] #added to prevent MaxIntIX does not exist errors
    IntBasin_labels = skimage.measure.label(IntBasin)
    for r in np.arange(1,np.max(IntBasin_labels)):
        PixelIdxList = np.argwhere(IntBasin_labels==r)
        ixm = np.argmax(DD[IntBasin_labels==r]);
        MaxIntIX = PixelIdxList[ixm];

        Iobj[PixelIdxList[0][0],PixelIdxList[0][1]] = 0;
        SILLS[PixelIdxList[0][0],PixelIdxList[0][1]] = 1;
    ixm = MaxIntIX;
    Iobj[ixm[0],ixm[1]] = 0;
    SILLS[ixm[0],ixm[1]] = 1;

    # establish the connectivity between sills and flats
    #dem = ZintoDB;
    whereSILLS = np.argwhere(SILLS);
    rows=[]
    cols=[]
    for rowcol in whereSILLS:    
        [row,col] = rowcol
        rows = np.append(rows,row)
        cols = np.append(cols,col)

    IXsill    = [rows,cols];
    rowadd = [-1, -1, 0, 1, 1,  1,  0, -1];
    coladd = [ 0,  1, 1, 1, 0, -1, -1, -1];
    PreSillPixel = [0]
    for r  in np.arange(0,8):
        rowp = rows + rowadd[r];
        colp = cols + coladd[r];

        ValidRowColPair1 = np.bitwise_and(rowp>0, colp>0)
        ValidRowColPair2 = np.bitwise_and(rowp<Z.shape[0], colp<Z.shape[1])
        ValidRowColPair  = np.bitwise_and(ValidRowColPair1, ValidRowColPair2) 
        whereValidRowColPair = np.where(ValidRowColPair)

        IXPreSill = [rowp[whereValidRowColPair],colp[whereValidRowColPair]];
        I1 = np.ravel_multi_index([np.int_(rows[whereValidRowColPair]), np.int_(cols[whereValidRowColPair])],Z.shape)
        I2 = np.ravel_multi_index([np.int_(IXPreSill[0]),np.int_(IXPreSill[1])],Z.shape)
        I3 = np.ravel_multi_index([np.int_(IXPreSill[0]),np.int_(IXPreSill[1])], Z.shape)
        #PreSillPixelCondition = (np.argwhere(np.bitwise_and((Z_ravel[I1] == 
         #           Z_ravel[I2]),
         #          Iobj.ravel()[I3]))         
         #   if np.count_nonzero(PreSillPixelCondition)>0:
         #       for i in  np.arange(0,len(PreSillPixelCondition)):
         #           PreSillPixelAddition = ([IXPreSill[0][PreSillPixelCondition[i]],IXPreSill[1][PreSillPixelCondition[i]]])
         #           PreSillPixel.append(PreSillPixelAddition)
         #   else:
         #       continue
    PreSillPixel.pop(0);
        

    Iobj  = np.bitwise_not(Iobj)    
    D = scipy.ndimage.distance_transform_edt(Iobj)
    masked = np.inf * np.ones(Z.shape,D.dtype); 
    masked[Iobj] = 0;
    D[np.bitwise_not(Iobj)]=np.inf
    D = ((skimage.morphology.reconstruction(seed = D+1,mask=masked,method='erosion'))- D) *cellsize
    D = np.nan_to_num(D)   

    D[Iobj] = 0
    #D = D**-1
    cost_field  = map_image_to_costs(D**-1,PreSillPixel)
    D = _wdt_python(cost_field) +1
    D[Iobj] = -np.inf
    
    #del PreSillPixel
    V = np.reshape(D.data,[1,D.shape[0]*D.shape[1]])
    if np.any(np.isnan(np.diff(D.ravel()))):
        IXSortedFlats = np.arange(0,len(Z_ravel))
        IXSortedFlats = IXSortedFlats[::-1]
    else:
        IXSortedFlats = np.argsort(D.ravel());
        IXSortedFlats = IXSortedFlats[::-1]
    #del D

    ndx = np.arange(np.uint32(0),np.uint32(nrc));
    ndx = ndx[IXSortedFlats];
    
    ndx = np.arange(np.uint32(0),np.uint32(nrc));
    ndx = ndx[IXSortedFlats];
    del IXSortedFlats

    ix = np.argsort(Z_ravel[ndx]);
    ix = ix[::-1]
    ix = ndx[ix]
    del ndx

     # a fast solution that has quite much memory overhead...
    pp = np.zeros(Z_ravel.shape,dtype=np.int32);
    IX = np.arange(np.int32(0),np.int32(Z_ravel.shape));
    pp[ix] = IX;
    pp = pp.reshape(Z.shape)
    
    # cardinal neighbors
    IXC1 = skimage.morphology.dilation(pp, skimage.morphology.selem.diamond(1))
    IXC1 = IXC1.ravel()
    xxx1 = IXC1;
    IX   = IXC1[ix];
    IXC1 = ix[IX];
    
    G1   = (Z_ravel[ix]-Z_ravel[IXC1])/(cellsize);
    I4 = (np.argwhere(ix == IXC1)).ravel()

    I4 = list(I4)
    I4_test = np.zeros(G1.shape)
    I4_test[I4] = -np.inf
    G1 = G1 + I4_test;
    #G1[ix == IXC1] = -np.inf;

     # diagonal neighbors
    kernel = np.array([[1,0,1],[0,1,0],[1,0,1]])
    IXC2 = skimage.morphology.dilation(pp,kernel);
    IXC2 = IXC2.ravel()
    xxx2 = IXC2;
    IX   = IXC2[ix];
    IXC2 = ix[IX];
    G2   = (Z_ravel[ix]-Z_ravel[IXC2])/np.linalg.norm([cellsize,cellsize]);


    # choose the steeper one
    #I  = np.bitwise_and(G1<=G2, xxx2[ix]>xxx1[ix]);
    I  = dask.array.bitwise_and(dask.array.less_equal(G1,G2),xxx2[ix]>xxx1[ix])
    ixc = IXC1;
    ixc[I] = IXC2[I];

    I = ixc == ix;
    ix = ix[np.bitwise_not(I)];
    ixc = ixc[np.bitwise_not(I)];

    # remove nans
    I = np.isnan(Z_ravel);
    ixc = ixc[~I[ix]];
    ix = ix[~I[ix]];
    
    ix = np.int_(ix[~np.isnan(ix)])
    ixc = np.int_(ixc[~np.isnan(ixc)])
    
    DBcounter = 0;
    D = np.zeros(Z_ravel.shape[0],dtype=np.int32);
    outlets=np.zeros((len(ix),1))
    for r in np.arange(len(ix)-1,1,-1):
        if D[ixc[r]] == 0:
            DBcounter = DBcounter+1;
            D[ixc[r]] = DBcounter;
            outlets[DBcounter] = ixc[r];

        D[ix[r]] = D[ixc[r]];

    D = D.reshape(Z.shape)
    return D

def cleandrainagebasins(D):
    D_labelled = skimage.morphology.label(D)

    D_labelled_new = D

    #merge down
    for i in np.arange(25,D.shape[0]-1,25):
        topside = D[i-1,:]
        topside_values = np.argwhere(np.diff(topside)!=0)
        for ii in np.arange(0,len(topside_values)):
            col_to_merge = topside_values[ii]
            basin_to_merge = D_labelled[i,col_to_merge]
            mask_right = np.argwhere(D_labelled == basin_to_merge)
            [mask_x,mask_y] = np.array(mask_right).T
            D_labelled_new[mask_x,mask_y] = D[i-1,topside_values[ii]]

    #Merge left
    for i in np.arange(25,D.shape[1],25):
        leftside = D_labelled[:,i-1]
        leftside_values = np.argwhere(np.diff(leftside)!=0)
        for ii in np.arange(0,len(leftside_values)):
            row_to_merge = leftside_values[ii]
            basin_to_merge = D_labelled[row_to_merge,i]
            D_labelled_new[D_labelled == basin_to_merge] = D[leftside_values[ii],i-1]

    xarray_D = xr.DataArray(data = D_labelled_new, coords = sample.coords,dims = sample.dims, attrs = sample.attrs)
    return xarray_D

def Inan(dem):
    Inan = np.isnan(dem)
    return Inan

def filledbasins(dem,Inans):
    marker = np.negative(dem);
    #dem = dem + np.multiply(Inans,-np.inf)
    II = np.zeros(dem.shape);
    II[1:-1,1:-1] = 1;
    mask = np.int_(np.bitwise_and(np.bool_(II),np.bitwise_not(np.bool_(Inans))))*-np.inf
    marker = marker+mask;
    demfs = -skimage.morphology.reconstruction(marker,-dem, method='dilation')
    return demfs