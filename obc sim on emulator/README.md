This folder contains files regarding running the simulation on obc. 

`dummy_simobc_data.py` was made to generate representable data in form of format to test the parser and plotter. 

`obc_em_parser.py` is the parser that goes through data collected on putty and collects data relevant to us in seperate txt files. To use it, you can just run it and the default input file is out.txt and default output files are parsed_q/ohm/com/eclipse.txt. If you want them to be something else, you can run this as 
```
python .\obc_em_parser.py inputfile.txt outputfile
``` 

Note there is no .txt in the output file.  

`plot_obc_em.py` plots the results from the data collected in the txt files. One can just run it and the default input files for this are parsed_q/ohm/com/eclipse.txt. If you want to plot data that is saved in different txt files, you'd run this as 
```
python .\plot_obc_em.py inputfile
```

Again, note no .txt in the name.

The folder also contains some parsed data from one of the collected logs and some example plots. 
