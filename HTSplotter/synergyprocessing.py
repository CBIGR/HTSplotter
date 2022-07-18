import numpy as np
import os
import math

class Commonprocessing():
    def __init__(self, grup):
        self.grup = grup
        self.name_list = []

    def get_comb_name(self):
        counter = None

        for i in range(self.grup[0].concentration.shape[0]):
            fr = ''
            fr, counter = self.name_rec(fr, self.grup[1:], counter)

            self.name_list.append(fr)

        return self.name_list

    def name_rec(self, fr, grup, counter, ind=0):

        if len(grup) > ind:
            if counter is None:
                counter = np.zeros(len(grup), dtype=np.int16)

            group = grup[ind]

            for i in range(ind, -1, -1):
                if counter[i] == grup[i].concentration.shape[0]:
                    counter[i] = 0
                    if i > 0:
                        counter[i - 1] += 1

            fr, counter = self.name_rec(fr, grup, counter, ind + 1)

            fr = str(group.name) + " " + str(group.concentration[counter[ind], 0]) + ' ' + str(group.unidade) + '_' + fr

            if ind == len(counter) - 1:
                counter[ind] += 1

            return fr, counter

        else:

            return fr, counter

    @staticmethod
    def adjust_data_new(data):

        for i in range(len(data)):
            for j in range(len(data[i].data_grup)):
                for k in range(len(data[i].data_grup[j])):
                    if data[i].data_grup[j, k] < 0:
                        data[i].data_grup[j, k] = 0

    @staticmethod
    def adjust_combinationvalues(data_grup):
        for i in range(len(data_grup[0].data_grup)):
            for k in range(len(data_grup[0].data_grup[i])):
                if data_grup[0].data_grup[i, k] < 0:
                    data_grup[0].data_grup[i, k] = 0