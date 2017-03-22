# Motif similarity

We estimate the similarity and its significance (Q-value) between two sets of motifs. The similarity is taken to be one minus the Shannon-Jensen distance. More specifically we follow the procedure presented in Itzkovitz et al. 2006, "Coding limits on the number of transcription factors".

Example: 

PWM_motif_set1 =  {  
[0.17 0.00 0.56 0.17 0.23 0.00 0.00 0.00  
 0.83 0.10 0.35 0.00 0.00 0.15 1.00 0.08  
 0.00 0.55 0.09 0.83 0.08 0.15 0.00 0.92  
 0.00 0.35 0.00 0.00 0.69 0.70 0.00 0.00],  
[0.00 0.00 0.74 0.67 0.00 0.00 0.08 0.00  
 1.00 0.00 0.00 0.00 0.34 0.28 0.36 0.00  
 0.00 1.00 0.26 0.10 0.53 0.00 0.56 1.00  
 0.00 0.00 0.00 0.23 0.13 0.72 0.00 0.00]  
};

PWM_motif_set2 =  {  
[0.17 0.00 0.56 0.17 0.23 0.00 0.00 0.00  
 0.83 0.10 0.35 0.00 0.00 0.15 1.00 0.08  
 0.00 0.55 0.09 0.83 0.08 0.15 0.00 0.92  
 0.00 0.35 0.00 0.00 0.69 0.70 0.00 0.00],  
[0.00 0.00 0.74 0.67 0.00 0.00 0.08 1.00  
 0.00 0.50 0.00 0.00 0.34 0.28 0.36 0.00  
 0.50 0.00 0.26 0.10 0.53 0.00 0.56 0.00  
 0.50 0.50 0.00 0.23 0.13 0.72 0.00 0.00]  
};

[\~,\~,\~,qval] = measure_motifs_similarity(PWM_motif_set1,PWM_motif_set2)

qval =

    0.0003    0.0006
    0.0005    0.0057
