###! CALC LD DECAY??
# https://eacooper400.github.io/gen8900/exercises/ld-2.html


# FROM THE WZA PAPER:
# LD decayed rapidly in our simulations with SNPs that were ~ 600 bp apart having, on average, half the LD of immediately adjacent SNPs in neutral simulations




### Find the point when LD drops below 0.1
### The which statement tells me which rows of the table have means
### greater than 0.1,
### the brackets are giving me the subset of the "bins" with averages
### that correspond to these rows, and the max statement tells me which bin has
### the greatest distance value
LD.drop = max(LD.averages$bins[which(LD.averages$my.means>0.1)])
abline(v=LD.drop, col="blue", lty=2) # Add a vertical line corresponding to the drop off point
text(25500, 0.8, "LD drops below 0.1", col="blue") # label that line




### Find the LD half life
LD.half = (max(LD.averages$my.means))/2 # Calculate half of the maximum average LD value
LD.half.point = max(LD.averages$bins[which(LD.averages$my.means>LD.half)]) # Same strategy as when I found the 0.1 drop off point
abline(v=LD.half.point, col="darkgreen", lty=2) # add a vertical line at the half life point
text(8500, 0.9, "LD half life", col="darkgreen") # label this line

