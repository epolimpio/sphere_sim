from myutils import extractVariablesFromString

string = 'abcd_04-bb_kc-99-00'
pattern = 'abcd_{0}-{1}_{2}-{3}-{4}'

extractVariablesFromString(string,pattern)