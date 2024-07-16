# Maximum Defective Clique

Maximum defective clique algorithm 

## Build
```bash
cmake . && make
```


## Run
- Command
```bash
bin/run -d <dataset> -k <k> -a <algorithm>

usage: bin/run -d=string -k=int [options] ...
options:
  -d, --data    dataset path (string)
  -k, --key     value of k (int)
  -a, --algo    algorithm (string [=MDC])
  -h, --help    print this message
```
> Available algoritms: MDC / RussianDoll / KDBB

- Example
```bash
bin/run -d datas/socfb-Harvard1 -k 1 -a MDC
```

- Data Format

The input data should be given as a list of edges, which follows the format below: 
```
<number of edges> <number of L vertices> <number of R vertices>
<v1> <v2>
<v3> <v4>
...
```
