## Running ceif
ceif is a command-line program controlled by arguments, basic syntax is:

    ceif [OPTION]...

Input data is assumed to be comma separated values. Different separator can be given by option -f.
### Options

| Option | Purpose  |
|:----|----|
| -h&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| Display help and exit|
| -d&nbsp;INTEGER | Number of decimals when printing and saving dimension values. Default is 6|
| -V | Display version and exit|
| -I&nbsp;LIST | LIST is a comma separated list field numbers (first = 1) which should be ignored when reading the input file. Ranges can be given with dash (e.g. 2-9). Default is to read all fields|
| -U&nbsp;LIST | LIST is a comma separated list field numbers (first = 1) which should be processed when reading the input file. Ranges can be given with dash (e.g. 2-9). This can be used with -I when there are lot of unused fields|
| -t&nbsp;INTEGER&nbsp;&nbsp; | Number of trees to use. Default is 100|
| -s&nbsp;INTEGER | Number of samples for each tree. Default is 256|
| -f&nbsp;CHAR | Field separator for input files|
| -l&nbsp;FILE | File to be used in algorithm training| 
| -a&nbsp;FILE | File to analyse|
| -c&nbsp;FILE | File to categorize|
| -p&nbsp;STRING | Printf style format to print anomaly data or categorized data. See printing directives below|
| -o&nbsp;FILE | Print output to FILE. Default is to use stdout|
| -w&nbsp;FILE | Write forest data to FILE. Typically result of analysing data using file with option -l. Data can be later read with option -r|
| -O&nbsp;FLOAT| Outlier score for anomaly detection. Data with higher or equal score is considered as an anomaly and printed with format given by option -p. Use values 0.0 - 1.0|
| -r&nbsp;FILE | Read forest data from file. File should have been written earlier with option -w|
| -C&nbsp;LIST | List of field numbers to be used as a category field. Default is not to use category field. Field values are separated by colon to form a category string|
| -L&nbsp;LIST | List of field numbers to be used as a label field. Default is not to use label field. Field values are separated by colon to form a label string|
| -F&nbsp;REGEXP | Filter categories using regular expression. Forests having category string matching REGEXP are not used in analysis or categorization. Several options can be given. If REGEXP is preceded by "-v " then matching is inverted|
| -H | Input data contains a header line which is ignored. Default is to read all lines|
| -S | Set locale to local locale. Default is use locale "C"|
| -R&nbsp;FLOAT| Interception point ***p*** range extension factor. Higher value means that starting interception points are selected more away from sample data points and maximum tree height is larger|
| -T&nbsp;FLOAT| Generate test data. Test data is generated using sample set min/max values and test data point interval given by option -i (default is 256). Test data range can be enlarged by FLOAT. E.g. value 1.0 doubles the test data range. After test data is printed max 10240 sample data points are printed with score 0|
| -i&nbsp;INTEGER| Test data point interval. Larger value means more dense test data point set|
| -u&nbsp;INTEGER| Accept only unique samples when sampling input data. INTEGER is value between 0..100 (default is 10). This is the percentage of input data rows to be checked for uniqueness. Value 100 can be used if every accepted sample data should be unique|
| -m&nbsp;STRING| Printf format for printing float values for sample and sample average values. Default is "%.*f"|
| -e&nbsp;CHAR| Value separator when printing sample, sample average and analysed data values. Default is comma|
| -n | adjust n-vector to be perpendicular to dimension attribute having largest value range|
| -M&nbsp;STRING | Print category value, average values or last update time of forests which have not used in analysis. Optional printf format STRING is used in printing|
| -D&nbsp;INTEGER | Before saving the forest data to file delete forests which have not been updated INTEGER (seconds) ago. If INTEGER is followed by a letter from set Y,M,D or m the INTEGER is consired to be years, months, days or minutes.|

If FILE is "-" then standard input or output is read or written.

#### Printing directives

| Directive | Meaning |
|----|----|
| %r | Current input file row number|
| %s | Anomaly score |
| %c | Category values separated by colon from input data|
| %l | Label value from input data|
| %d | Separated list of dimension values|
| %a | Separated list of dimension average values|
| %v | Current input row values|
| %x | Outlier score value in RGB values. Presented as hex value (.e.g 127F77). Note that value zero is printed as black|
| %C | Found category values separated by colon when categorizing data|
| %t | Time when category has been last updated. In human readable form using current locale|
| %: | Category value separator|
| %. | Label value separator|
| %% | Percent sign|

Value separator for d,a and v can be given by option -e.

### Examples

#### Learn and analyze same field
Learn and analyse file data.csv, use fields 1-3 (by ignoring fields 4-100) for analysis and print anomalies having score 0.6 or greater.

    ceif -l data.csv -a data.csv -R30 -I4-100  -O0.6

#### Learn and write forest data to file data.f

    ceif -l data.csv -w data.f -R30 -I4-100 -O0.6

#### Read forest data from data.f and analyse data.csv
Note that parameters like -R and -O are saved to forest file and need not be given again.  

    ceif -r data.f -a data.csv

####  Read forest data from data.f and add more samples from data2.csv 
Note that -w must given in order to save enhanced forest data.

    ceif -r data.f -l data2.csv -w data.f

#### Learn and write forest data with category
Use field number 5 as category field

    ceif -l data.csv -w data.f -R30 -I4-100 -O0.6 -C 5

#### categorize data from data.csv using forest data from previous example
Print analyzed category value (%C) and select fields as comma separated list (%v) for each input row from data.csv.

    ceif -r data.f -c data.csv -p "%C %d"

#### Generate test data set using forest data file
Generate data set around sample data points by enlarging the area with factor 1. Every dimension value range consists of 512 test values. 
Test data sample values and outlier score value in RGB value separated by semicolon are printed to file plot_data.csv.

    ceif -r data.f -T1 -i512 -e";" -p"%d;0x%x" -o plot_data.csv
