#!/usr/bin/python3

file = open('tree.nwk', 'r')
out = open('yuyt2.txt', 'w')

line = file.readline()
file.close()
#print(line.count('\'')/2)

line = line.replace('\'', '')
line = line.replace('(', '')
line = line.replace(')', '')
line = line.replace(' ', '')
proteins = line.split(',')
for p in proteins:
    #out.write('_'.join(p.replace('_', ' ').split(' ')[0:2]) + '\n')
    out.write('_'.join(p.split('_')[0:2]) + '\n')
out.close()
