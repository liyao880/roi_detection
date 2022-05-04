__author__ = "Sherif Farag"
__email__ = "sherif_farag@med.unc.edu"
__status__ = "Development"

import argparse
import numpy as np
from PIL import Image
import h5py
import os


def normalizeStaining(img, Io=240, alpha=1, beta=0.15):
    ''' Normalize staining appearence of H&E stained images
    
    Example use:
        see test.py
        
    Input:
        I: RGB input image
        Io: (optional) transmitted light intensity
        
    Output:
        Inorm: normalized image
        H: hematoxylin image
        E: eosin image
    
    Reference: 
        A method for normalizing histology slides for quantitative analysis. M.
        Macenko et al., ISBI 2009
    '''
    
    #print("noussa: ", img.shape)
    HERef = np.array([[0.5626, 0.2159],
                      [0.7201, 0.8012],
                      [0.4062, 0.5581]])
        
    maxCRef = np.array([1.9705, 1.0308])
    
    # define height and width of image
    h, w, c = img.shape
    #print("height: ",h,"width: ",w,"Channel: ",c)
    
    # reshape image
    img = img.reshape((-1,3))

    # calculate optical density
    OD = -np.log((img.astype(np.float)+1)/Io)
    
    # remove transparent pixels
    ODhat = OD[~np.any(OD<beta, axis=1)]
        
    # compute eigenvectors
    eigvals, eigvecs = np.linalg.eigh(np.cov(ODhat.T))
    
    #eigvecs *= -1
    
    #project on the plane spanned by the eigenvectors corresponding to the two 
    # largest eigenvalues    
    That = ODhat.dot(eigvecs[:,1:3])
    
    phi = np.arctan2(That[:,1],That[:,0])
    
    minPhi = np.percentile(phi, alpha)
    maxPhi = np.percentile(phi, 100-alpha)
    
    vMin = eigvecs[:,1:3].dot(np.array([(np.cos(minPhi), np.sin(minPhi))]).T)
    vMax = eigvecs[:,1:3].dot(np.array([(np.cos(maxPhi), np.sin(maxPhi))]).T)
    
    # a heuristic to make the vector corresponding to hematoxylin first and the 
    # one corresponding to eosin second
    if vMin[0] > vMax[0]:
        HE = np.array((vMin[:,0], vMax[:,0])).T
    else:
        HE = np.array((vMax[:,0], vMin[:,0])).T
    
    # rows correspond to channels (RGB), columns to OD values
    Y = np.reshape(OD, (-1, 3)).T
    
    # determine concentrations of the individual stains
    C = np.linalg.lstsq(HE,Y, rcond=None)[0]
    
    # normalize stain concentrations
    maxC = np.array([np.percentile(C[0,:], 99), np.percentile(C[1,:],99)])
    tmp = np.divide(maxC,maxCRef)
    C2 = np.divide(C,tmp[:, np.newaxis])
    
    # recreate the image using reference mixing matrix
    Inorm = np.multiply(Io, np.exp(-HERef.dot(C2)))
    Inorm[Inorm>255] = 254
    Inorm = np.reshape(Inorm.T, (h, w, 3)).astype(np.uint8)  
    
    # unmix hematoxylin and eosin
    H = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,0], axis=1).dot(np.expand_dims(C2[0,:], axis=0))))
    H[H>255] = 254
    H = np.reshape(H.T, (h, w, 3)).astype(np.uint8)
    
    E = np.multiply(Io, np.exp(np.expand_dims(-HERef[:,1], axis=1).dot(np.expand_dims(C2[1,:], axis=0))))
    E[E>255] = 254
    E = np.reshape(E.T, (h, w, 3)).astype(np.uint8)
    
    # if saveFile is not None:
    #     Image.fromarray(Inorm).save(saveFile + "/" + str(patch_number) + '.png')
        #Image.fromarray(H).save(saveFile+'_H.png')
        #Image.fromarray(E).save(saveFile+'_E.png')

    return Inorm.astype(np.uint8), H, E
    

def computeH5(f, save_name):
    file = h5py.File(f, 'r')
    dset = file['imgs']
    coords = file['coords']
    print("Number of patches:")
    print("dset.shape: ", dset.shape, "coords.shape: ", coords.shape)
    # if f.startswith("_"):
    #     n=f.split("_")[2]
    #     pos=int(f.split("_")[3].split(".h5")[0])
    # else:
    #     n=f.split("_")[-2]
    #     pos=int(f.split("_")[-1].split(".h5")[0])
    #n=len([name for name in os.listdir(args.saveFile) if os.path.isfile(args.saveFile + "/" + name)])
    # print("n: ", n)
    # print("pos: ", pos)
    error_l=[]
    img_norm = []
    coords_norm = []
    for i in range(0, len(dset)):
        if(coords[i][0]==34752 and coords[i][1]==25344):print("eshta")
        img=dset[i]
        try:
            #print(i, coords[i])
            Inorm, _, _ = normalizeStaining(img = img, Io = args.Io, alpha = args.alpha,beta = args.beta)
            img_norm.append(Inorm)
            coords_norm.append(coords[i])
        except:
            print(i, coords[i])
            error_l.append(i)
        # print("#####################")
        # print()
    img_norm = np.stack(img_norm)
    coords_norm = np.stack(coords_norm)
    print("douds: ", len(error_l))
    if np.all(np.isfinite(img_norm)):
        save_file = h5py.File(save_name, "w")
        aset = save_file.create_dataset('imgs', shape=img_norm.shape, dtype='int8')
        aset[:] = img_norm
        aset.attrs['level_dim'] = dset.attrs['level_dim']
        file.close()
        cset = save_file.create_dataset('coords', shape=coords_norm.shape, dtype='int32')
        cset[:] = coords_norm
        save_file.close()
    else:
        file.close()
    
    
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--imageFile', type=str, default=r'D:\yaoli\data\melanoma\patches_all\patches\PGEM1179.h5', 
                        help='/pine/scr/y/a/yaoli/data/melanoma/patches_all/patches_norm/splitPatches/PGEM862/PGEM862.h5_5000_0.h5, RGB image file')
    parser.add_argument('-o', '--saveFile', type=str, default=r'D:\yaoli\data\melanoma\patches_all_norm\patches\PGEM1179.h5', help='/pine/scr/y/a/yaoli/data/melanoma/patches_all/patches_norm/splitPatches/PGEM862/results/norm, save file')
    parser.add_argument('--Io', type=int, default=240)
    parser.add_argument('--alpha', type=float, default=1)
    parser.add_argument('--beta', type=float, default=0.15)
    args = parser.parse_args()
    computeH5(args.imageFile, args.saveFile)
    print("DONE!")

    
