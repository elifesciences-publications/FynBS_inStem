'Directory'[Y:\Imag_Data\U2OS_Cell\] is the path where list of cells in separate folder[Cell(Number), say Cell07] is saved with the
Donor Acceptor and FRET_ch images, names containing the channel name. The image format is .tif

0. Copy the FRET_ch image as a mock file and rename it as FRET_Index.[This is just a mock file]
1. Run calc_fret_index. moth_dirt should be the directory path. 
Input of folder name[fld_nm] will be asked. Press Enter followingly.
2. From here onwards moth_dirt variable should be 'Directory\fld_nm'
3. Run cell_intensity_quantification. Generates the Quandrant Data for all channels.

4. Put a .mat file named bef_idx.mat containing a row array before_idx in the path 'Directory\fld_nm' , which stores
   sequentially the point of PDGF addition for following computations

5. Run HQ_fret_parameters. High Quadrant - parameters
6. Run initial_int_quant. Initail Intensity of the whole cell.
7.1 Run Line_Scan. Saves the line scan data
7.2 Run Line_Scan_Analysis. Performs the final quantification and saves the quantified data
8.1 Run morphology_changes. Performs the data quantification under sub folder Morph_Quant
8.2 Run morph_fret_corr. Input a 0 or 1 based on the requirement to perform the computation on Before or
After.
Performs the corealtion under sub folder Morph_Quant
9. Run polr_extoscc_pulst. The final quantification.