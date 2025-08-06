# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 00:11:47 2023

@author: Yvan Guo

"""

import math

def purity(a, b):
    a = set(a)
    b = set(b)
    # calucate jaccard similarity
    j = float(len(a.intersection(b)))/len(a)
    return j

def jaccard_similarity(a, b):
    # convert to set
    a = set(a)
    b = set(b)
    # calucate jaccard similarity
    j = float(len(a.intersection(b))) / len(a.union(b))
    return j


def GMeasure(a, b):
    # convert to set
    a = set(a)
    b = set(b)
    # calucate jaccard similarity
    j = float(math.sqrt(len(a.intersection(b))/len(a)*len(a.intersection(b))/len(b)))
    return j


def main():
      ...

if __name__ == '__main__':
    main()