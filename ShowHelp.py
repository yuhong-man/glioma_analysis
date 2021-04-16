#!/data2/wangb/bin/python3
'''
Show the content after the last __main__. 
'''
import sys 

def printAfterLastMain(file_in, keyword=''):
    print('####################Help:')
    main_passed = False
    for line in open(file_in):
        line = line.strip()
        if '__main__' in line:
            main_passed = True 
        if main_passed and keyword in line:
            print(line)
    print('####################HelpEnd')

printAfterLastMain(sys.argv[0])

        

