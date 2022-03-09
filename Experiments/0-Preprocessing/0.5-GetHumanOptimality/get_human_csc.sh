#!/bin/bash

# script to collect the codon stability coefficient data from the paper: 
# "Translation affects mRNA stability in a codon-dependent manner in human cells"
# produced by the Bazzini lab 

# make dir to put data 
mkdir -p data

# download the dataset 
wget https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNDUzOTYvZWxpZmUtNDUzOTYtZmlnMS1kYXRhMi12Mi5jc3Y-/elife-45396-fig1-data2-v2.csv?_hash=oV0Fjo95uQOzu5LreFXU9sbiAG2ub8ZzXLyP%2B0iTk98%3D
mv elife-45396-fig1-data2-v2.csv\?_hash\=oV0Fjo95uQOzu5LreFXU9sbiAG2ub8ZzXLyP+0iTk98\= data/human_csc_all_methods.csv

