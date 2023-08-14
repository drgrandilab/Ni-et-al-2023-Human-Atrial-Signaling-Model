### Matlab script to analyze AP and Ca instability, specifically DADs and SCR in each cell.


- To see if the code works, modify line 16 of file process_folder.m

``` matlab
for id = 0:599
```
to 

```matlab 
for id = 0:1
```

and run 

```matlab 
process_folder('./') 
```
in Matlab command window