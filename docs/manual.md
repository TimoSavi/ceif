## Running ceif
ceif is a command-line program controlled by arguments, basic syntax is:

    ceif [OPTION]...

Input data is assumed to be comma separated values. Different separator can be given by option -f.
### Options

| Option | Purpose  |
|:----|----|
| -h&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| Display help and exit|
| -d&nbsp;INTEGER | Number of decimals when printing and saving dimension values. Default is 6|
| -V | Display version and exit|
| -I&nbsp;LIST | LIST is a comma separated list of field numbers (first = 1) which should be ignored when reading the input file. Ranges can be given with dash (e.g. 2-9). Default is to read all fields|
| -U&nbsp;LIST | LIST is a comma separated list of field numbers (first = 1) which should be processed when reading the input file. Ranges can be given with dash (e.g. 2-9). This overrides equal values from option -I|
| -X&nbsp;LIST | LIST is a comma separated list of field numbers (first = 1) which should be as text fields when reading the input file. Ranges can be given with dash (e.g. 2-9). A hash value from range 0-32770 is generated using field string value. Note that this is not collision free. This should be used mainly for simple classifications like "yes/no" or "Male/Female/Unknown". Note also that "Yes" and "yes" will produce different value|
| -G&nbsp;LIST | LIST is a comma separated list of dimensions attribute indices. Combination of these dimension attributes must have outlier score along total score. Ranges can be given using dash. These are dimensions attribute indices, not input line indices (first = 1)|
| -t&nbsp;INTEGER&nbsp;&nbsp; | Number of trees to use. Default is 100|
| -s&nbsp;INTEGER | Number of samples for each tree. Default is 256|
| -f&nbsp;CHAR | Field separator for input files|
| -l&nbsp;FILE | File to be used in algorithm training| 
| -a&nbsp;FILE | File to analyse|
| -c&nbsp;FILE | File to categorize|
| -p&nbsp;STRING | Printf style format to print anomaly data or categorized data. See printing directives below|
| -o&nbsp;FILE | Print output to FILE. Default is to use stdout|
| -r&nbsp;FILE | Read forest data from file. File should have been written earlier with option -w|
| -w&nbsp;FILE | Write forest data to FILE. Typically result of analysing data using file with option -l. Data can be later read with option -r|
| -z&nbsp;FILE | Read and write forest data from/to file. Forest data is read from file FILE and after any processing written back to FILE|
| -O&nbsp;FLOAT| Outlier score for anomaly detection. Data with higher or equal score is considered as an anomaly and printed with format given by option -p. Use values 0.0 - 1.0|
| -O&nbsp;FLOATs| Scaled outlier score for anomaly detection. The actual analyzed score is scaled to range 0-1 using forest min/max scores. This ensures that the best inlier has value zero and the farthest outlier will get value near 1.0. When given with categorize option (-c) the results are limited by this value. Only category results having lower value than this are accepted.  Use values 0.0s - 1.0s|
| -O&nbsp;FLOAT%| Outlier score for anomaly detection is calculated using sample score distribution taking the score value under which FLOAT percent of samples have lower score. Use values  0% - 100%|
| -C&nbsp;LIST | List of field numbers to be used as a category field. Default is not to use category field. Field values are separated by colon to form a category string|
| -L&nbsp;LIST | List of field numbers to be used as a label field. Default is not to use label field. Field values are separated by colon to form a label string|
| -F&nbsp;REGEXP | Filter categories using regular expression. Forests having category string matching REGEXP are not used in analysis or categorization. Several options can be given. If REGEXP is preceded by "-v " then matching is inverted|
| -H | Input data contains a header line which is ignored. Default is to read all lines|
| -S | Set locale to local locale. Default is use locale "C"|
| -T&nbsp;FLOAT| Generate test data. Test data is generated using sample set min/max values and test data point interval given by option -i (default is 256). Test data range can be enlarged by FLOAT. E.g. value 1.0 doubles the test data range. After test data is printed max 10240 sample data points are printed with score 0|
| -i&nbsp;INTEGER| Test data point interval. Larger value means more dense test data point set|
| -u&nbsp;INTEGER| Accept only unique samples when sampling input data. INTEGER is value between 0..100 (default is 10). This is the percentage of input data rows to be checked for uniqueness. Value 100 can be used if every accepted sample data should be unique|
| -m&nbsp;STRING| Printf format for printing float values for sample and sample average values. Default is "%.*f"|
| -j&nbsp;STRING| Printf format for printing all single dimension attribute related metrics together. Metrics are printed using printf directive %m|
| -e&nbsp;CHAR| Value separator when printing sample, sample average and analysed data values. Default is comma|
| -M&nbsp;STRING | Print category value, average values or last update time of forests which have not used in analysis. Optional printf format STRING is used in printing|
| -D&nbsp;INTEGER | Before saving the forest data to file delete forests which have not been updated INTEGER (seconds) ago. If INTEGER is followed by a letter from set Y,M,D or m the INTEGER is consired to be years, months, days or minutes.|
| -N&nbsp;STRING | Print input values which are not assosiated with any categories. This can be used for printing "new" category values. Optional printf format STRING is used in printing|
| -A | Instead taking samples as they are, aggregate new samples values for each forest. Only one new aggregated sample for each forest is added for each usage of -l option|
| -q | Print forest information in human readable form and exit|
| -y | Print forest information ascii density map|
| -yy | Print forest information ascii density map with common sample scale for all forests|
| -E | Print samples with sample score|
| -k | Remove the sample having maximun sample score for each non filtered forest. If option is given several times, then several samples are removed. This can be used to remove outliers from samples. Modified sample set can be saved with option -w|
| -g&nbsp;FILE | Use the FILE as rc-file instead of ~/.ceifrc. Note that file given with option -g overrides options given before -g|
| -P | Print list of correlation coefficents with regression line slopea and y-intercepts for every dimension attribute pair and exit. Correlation coefficent is a value between -1.0 - 1.0|
| -R&nbsp;STRING | Remove all samples for a forest having STRING as forest string.|
| -v&nbsp;STRING | Print average score and other statistics calculated from analysed data using printing format STRING.|
| -Q&nbsp;STRING | Replace input data value using an expression in STRING, STRING is added to list of expression. If STRING starts with hyphen, then the expression is removed from the list. |


If FILE is "-" then standard input or output is read or written.

Default file format for options -r,-w and -z is JSON. If JSON is not available then CSV format is used. Ceif tries to obey the number of decimals (option -d) when saving data.
If no double formatting support is available, the number of decimals saved is the json library default.

#### Printing directives

| Directive | Meaning |
|----|----|
| %r | Current input file row number|
| %s | Anomaly score|
| %g | Anomaly score for dimensions given by option -G|
| %S | Average anomaly score for analysed data|
| %n | Number of rows for a forest|
| %o | Number of analyzed rows for a forest. This might be lower than %n value if data sampling is used (see ANALYZE\_SAMPLING in next section)|
| %h | Number of analysed rows having larger score than outlier score value|
| %c | Category string from input data. The original category when categorizing data|
| %C | Forest category string. The best matching category when categorizing data|
| %l | Label values|
| %d | Separated list of dimension values. Print text based dimensions as text|
| %u | Separated list of dimension values. Print text based dimensions as double|
| %a | Separated list of dimension average values|
| %e | Separated list of dimension attribute scores, ceif tries to analyze how each attribute affects the total score and gives each attribute a score|
| %m | Separated list of dimension metrics printed by attribute by attribute. Printf format for attribute metrics printed is given by option -j or rc-file variable PRINT\_DIMENSION. See rc-file for printing directives|
| %i | Dimension attribute index (first=1)|
| %v | Current input row values|
| %x | Outlier score value in RGB values. Presented as hex value (.e.g 127F77). Default gradient color for score values 0..1 is yellow to red. Colors can be changed using rc-file. Note that value zero is printed as black|
| %X | Outlier score value in RGB for dimensions given by option -G|
| %t | Time when category has been last updated. In human readable form using current locale|
| %: | Category value separator|
| %. | Label value separator|
| %% | Percent sign|

Value separator for d,u,a,m,e and v can be given by option -e.
Category value separator is semicolon and label value separator is dash. These can be changed in rc-file.

### User rc-file
Some defaults can be read from user specific rc-file ~/.ceifrc. File has variable-value pairs separated by whitespace. Comments start with #. 
These values can be overriden by command options and saved forest data (read using option -r). Note that file given with option -g overrides options given before -g.

Following variables are supported:

| Variable | Meaning | default value |
|----|----|----|
|SAMPLES|Number of samples taken for each tree, same affect as option -s|256|
|TREES|Number of trees for each forest, same affect as option -t|100|
|DECIMALS|Number of decimals used when saving forest data. Affects also printing of sample values (option -d)|6|
|AUTO\_SCALE|Scale sample values before analysing the forest, 1 = yes, 0 = no|1|
|CATEGORY\_SEPARATOR|Char to be used as a separator when concatenating category fields|;|
|LABEL\_SEPARATOR|Char to be used as a separator when concatenating label fields|-|
|OUTLIER\_SCORE|Outlier score for analysis, same values as for option -O can be used ("max", "average", float value 0..1 or float with suffix 's' 0s..1s)|
|MAX\_SAMPLES|Maximum number of samples for each forest|Default is calculated by number\_of\_trees * number\_of\_samples\_per\_tree|
|NEAREST|Score is adjusted by the distance to nearest sample point in leaf nodes, 1 = yes, 0 = no|1|
|ANALYZE\_SAMPLING|If the data to be analyzed is expected to be inpractical large it can be sampled. If the analyzed row count reaches the value defined by this variable then the sampling starts. Sampling is implement using reservoir sampling method. The number of analyzed rows is estimated to be k * (ln(x/k) + 1), where k = this parameter value, x = total row count|0 (default value, no sampling)|
|DEBUG|Print debug messages, 1 = yes, 0 = no|0|
|PRINT\_DIMENSION|Printf string for printf directive %m. This printf string can contain directives %d, %a, %e and %i ||
|DIM\_PRINT\_WIDTH|Attribute metrics printing width. Used when printing forest info with option -q. This can be used when number of dimension attributes is high and metrics do not fit to screen|25|
|CLUSTER\_SIZE|Ceif tries to find data cluster by taking the samples having lowest scores and counting the number of samples around them. Cluster size is fixed and is calculated by finding the distance from the sample having the lowest score to most distant sample. Cluster size is the distance multiplied by this value. Use values between 0-1|0.125|
|LOW_RGB_COLOR|RGB color code for score value 0 (printing directive %x). Values are given as hex string (e.g. 0x12fe44)|0xffff00, yellow|
|HIGH_RGB_COLOR|RGB color code for score value 1|0xff0000, red|

Example of rc-file:

    # My default values
    TREES 200
    CATEGORY_SEPARATOR +
    LABEL_SEPARATOR   .
    OUTLIER_SCORE 0.65s
    PRINT_DIMENSION "<value>%d</value>\n<score>%e</score>\n"

    #End of file

### Examples

#### Learn and analyze same file
Learn and analyse file data.csv, use fields 1-3 (by ignoring fields 4-100) for analysis and print anomalies having score 0.6 or greater.

    ceif -l data.csv -a data.csv -I4-100  -O0.6

#### Learn and write forest data to file data.f. Use scaled sample score value 0.5, 

    ceif -l data.csv -w data.f -I4-100 -O0.5s

#### Read forest data from data.f and analyse data.csv
Note that parameters like -O are saved to forest file and need not be given again.  

    ceif -r data.f -a data.csv

####  Read forest data from data.f and add more samples from data2.csv 
Note that -w must given in order to save enhanced forest data.

    ceif -r data.f -l data2.csv -w data.f

Option -z reads and saves to same file:

    ceif -z data.f -l data2.csv 

#### Learn and write categorized forest data 
Use field number 5 as category field

    ceif -l data.csv -w data.f -I4-100 -O0.6 -C 5

#### categorize data from data.csv using forest data from previous example
Print analyzed category value (%C) and select fields as comma separated list (%v) for each input row from data.csv.

    ceif -r data.f -c data.csv -p "%C %d"

#### Generate test data set using forest data file
Generate data set around sample data points by enlarging the area with factor 1. Every dimension value range consists of 512 test values. 
Test data sample values and outlier score value in RGB value separated by semicolon are printed to file plot\_data.csv.

    ceif -r data.f -T1 -i512 -e";" -p"%d;0x%x" -o plot_data.csv

#### Data value aggregation
If e.g. daily or hourly summaries should be analyzed then ceif should be called with -A option and each forest should have a daily or hourly aggregation key. 
Example file of hourly samples of internet traffic:

    $ cat traffic.csv
    Hour,Type,Inbytes,Outbutes
    12,In,123,1444
    12,Out,423,1644
    13,In,123,1444
    13,In,823,44
    14,In,13,1414
    14,In,9123,1444
    14,In,123,1443
    15,Out,12423,1644
    16,Out,423,16
    16,Out,493,1044
    17,Out,433,1644
    18,Out,493,1644

Making Hourly and traffic direction summaries, first two fields are category keys (-C1-2):

    ceif -l traffic.csv -C1-2 -H -A -d0 -w traffic.ceif

Forest file contents below, for each hour and traffic type the number of bytes are summarized:

    $ cat traffic.ceif
    G;2;"";"%s %v";100;256;"1-2";",";1;0.750000;10.000000;"";"";8;"";0;0;"";",";1;1
    F;"12:In";0.000000;0;1;1582453164
    S;123|1444
    F;"12:Out";0.000000;0;1;1582453164
    S;423|1644
    F;"13:In";0.000000;0;1;1582453164
    S;946|1488
    F;"14:In";0.000000;0;1;1582453164
    S;9259|4301
    F;"15:Out";0.000000;0;1;1582453164
    S;12423|1644
    F;"16:Out";0.000000;0;1;1582453164
    S;916|1060
    F;"17:Out";0.000000;0;1;1582453164
    S;433|1644
    F;"18:Out";0.000000;0;1;1582453164
    S;493|1644

And when the same is done next day and so on, a new sample is added for each forest (Hour and traffic direction combination).
After enough samples are collected a daily traffic can be analyzed for anomalies using -a option.

#### Using input data expressions (-Q)
Input data values can be modified before actual processing. Modification expression can have any tinyexpr expression (see https://github.com/codeplea/tinyexpr/tree/master?tab=readme-ov-file#grammar). Input fields are referred using '$n' notation where 'n' is the sequence number of the input data field (first field = 1). 

Expression has format:
```
$n=<expr>
```
or
```
$n=<expr>:d
```
```
'$n' is the reference to the field to be modified
'<expr>' is the tinyexp expression which can have '$n' references to input fields. Non existing fields are not replaced
':d' d is optional number of decimals to be used
```

Examples:

| Case  | Expression |
|:----|----|
| Divide the 4th field by 13, result has two decimals | '$4=$4/13:2' |
| Multiply the first field with second field | '$1=$1*$2' |
| Divide the second field with 60 and convert it to integer| '$2=floor($2/60)' |
| Set constant value for field 10 | '$10=3486' |

All expressions are saved in forest data. If an earlier saved expression should be removed, then prefix it with hyphen (e.g. -Q '-$4=$4/10') and save the data (-w or -z).

Example command line:
```
ceif -Q '$4=$4/13:2' -Q '$2=floor($2/60)' -Q '-$4=$4/10' ...
```
