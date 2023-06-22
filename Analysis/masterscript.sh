#!/bin/bash
myArray=(3_harmonisemultibi.R 7a_exploreMRresults.R 8a_createtables.R 8b_plot.R)
for str in ${myArray[@]}; do
chmod u+x ./$str
done
echo 'Initializing 3_harmonisemultibi.R' && ./3_harmonisemultibi.R && echo 'Initializing 7a_exploreMRresults.R' && ./7a_exploreMRresults.R && echo 'Initializing 8a_createtables.R' && ./8a_createtables.R && echo 'Initializing 8b_plot.R' && ./8b_plot.R && echo 'The master script finished without errors'
