## katherine carbeck
## 13 sept 2023
## run wza
## https://github.com/TBooker/WZA/blob/master/README.md

#columns in corr file:
# chrom  pos numPops maf TD_tau TD_pVal  MSP_tau MSP_pVal SHM_tau SHM_pVal DD_0_tau  DD_0_pVal  DD18_tau  DD18_pVal X bin_start bin_end window

# general_WZA_script.py 
    # --summary_stat = TD_pVal, MSP_pVal, SHM_pVal, DD_0_pVal, DD18_pVal


# TD
python general_WZA_script.py --correlations output_corr_windows.correlations.csv \
                              --summary_stat TD_pVal \
                              --window window \
                              --MAF maf \
                              --output TD_WZA.csv \
                              --sep "," 
# started running at 3:20 PM on 24 core machine 

# MSP
python general_WZA_script.py --correlations output_corr_windows.correlations.csv \
                              --summary_stat MSP_pVal \
                              --window window \
                              --MAF maf \
                              --output MSP_WZA.csv \
                              --sep "," &

# SHM_pVal                    
python general_WZA_script.py --correlations output_corr_windows.correlations.csv \
                              --summary_stat SHM_pVal \
                              --window window \
                              --MAF maf \
                              --output SHM_WZA.csv \
                              --sep "," &

# DD_0pVal
python general_WZA_script.py --correlations output_corr_windows.correlations.csv \
                              --summary_stat DD_0_pVal \
                              --window window \
                              --MAF maf \
                              --output DD_0_WZA.csv \
                              --sep ","  &                      

# DD18_pVal
python general_WZA_script.py --correlations output_corr_windows.correlations.csv \
                              --summary_stat DD18_pVal \
                              --window window \
                              --MAF maf \
                              --output DD18_WZA.csv \
                              --sep "," &
# started last 4 at 9:20PM

## _WZA.csv file columns:
# gene,SNPs,hits,Z,top_candidate_p,Z_pVal
    # gene - the name of the window
    # SNPs - the number of SNPs in this window
    # hits - the number of SNPs in the 99th percentile (not used for anything, just good to know)
    # Z - the Z score calculated for the gene
    # top_candidate_p - the result of the top-candidate method of Yeaman et al (2016 - Science)
    # LA - an indicator of whether the gene is causal for local adaptation
    # position - the average position of all SNPs in the window
    # Z_pVal - the p-value of the Z score (This is the WZA score)