import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import tifffile
from PIL import Image

from deepcell.utils.plot_utils import create_rgb_image, make_outline_overlay
from deepcell_toolbox import erode_edges

# download file
exec(open('https://github.com/KaWingLee9/in_house_tools/blob/main/multiplexed_images_pipeline/mesmer_function.py').read())

from deepcell.applications import Mesmer
app = Mesmer()

# The order of the image and the marker_list must be the same
marker_list = ['80ArAr','131Xe','SMA','Collagen_I','Vimentin','CD14','CK19','CD16','CD163','PanCK','CD15',
'PD-L1','CD31','CD45','CD57','CD117','FoxP3','CD4','E-cadherin','CD68','CD11c',
'CD20','CD8A','CD206','Arg1','PD-1','CD66b','GranzymeB','Ki67','CTLA-4','CD3',
'AFP','CD39','CD69','HLA-DR','GPC3','HepPar1','DNA1','DNA2','ICSK1','ICSK2','ICSK3']
# the keys of the dict can be whatever you like
in_nuc_marker_kwargs = {'nuc_channel': ['DNA1']}
in_mem_marker_kwargs = {# 'mem_channel_1': ['HepPar1','AFP','Arg1','GPC3'], 
                        # 'mem_channel_3':['SMA','Vimentin','CD31'],
                        # 'mem_channel_2': ['CD45','CD3', 'CD68',  'CD8A', 'CD66b'],
                        'mem_channel_1': ['ICSK1','ICSK2','ICSK3']
}

import os
tiff_path='/data/tiff/'
file_ls=os.listdir(tiff_path)

for file in file_ls:
    img_multipledxed=tifffile.imread(tiff_path+file)
    img_norm=Denoise_img(img_multipledxed,98)

    img_mesmer=Generate_mesmer_input(img_norm,marker_list,in_nuc_marker_kwargs,in_mem_marker_kwargs)
    # plot membrane channel and nuclear channel
    # Plot_mesmer_input(img_mesmer)

    # segmentation_predictions=app.predict(img_mesmer,preprocess_kwargs={'percentile':99.9},image_mpp=1,compartment='nuclear')
    try:
        segmentation_predictions=app.predict(img_mesmer,preprocess_kwargs={'percentile':99.9},image_mpp=1,compartment='whole-cell')
    except Exception:
        print(file)
        continue
    
    # plot membrane channel and nuclear channel with segmentation outline
    # rgb_img=create_rgb_image(img_mesmer,channel_colors=['blue', 'green'])
    # Plot_segmentation_result(rgb_img,segmentation_predictions)

    tifffile.imwrite('/data/'+file.split('.')[0]+'_mask.tiff',segmentation_predictions[0,:,:,0])