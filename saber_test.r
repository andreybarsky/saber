# test data filepaths, in different formats:

source("./saber_main.r")

text_data = 'test_data.txt' # assumes events separated by spaces, sequences separated by linebreaks, but you can specify alternatives for these
xl_data = 'test_data.xlsx' # assumes events separated by cells, sequences separated by rows
csv_data = 'test_data.csv'
whitespace_data = 'whitespace_data.txt'

txt_test = saber(text_data)
txt_test2 = saber(whitespace_data)
xl_test = saber(xl_data, min_obs=1) # sparse frequency matrix
csv_test = saber(csv_data, min_obs=1)
ws_test = saber(whitespace_data, lfe_collapse=0) # many low frequency events here