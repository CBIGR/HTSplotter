import numpy as np
import os
import math
from save_hdf5file import Hdf5database

class Blissmethod:

    def __init__(self, grup, data_grup, i, j):
        self.grup = grup
        self.data_grup = data_grup
        self.i = i
        self.j = j

        self.bivalues = []
        self.predvalues = []
        self.combname = []

        self.bliss_score()

    def bliss_score(self):

        pred_lis, bi_list, name_list = self.bi_score_matrix_new(self.grup, self.data_grup)
        self.predvalues = np.asarray(pred_lis)
        self.bivalues=np.asarray(bi_list)
        self.combname= name_list

    def bi_score_matrix_new(self, grup, data_grup):

        bi_list = []
        name_list = []
        # correct inhibition below control
        self.adjust_data_new(data_grup[1:])

        name_list = self.get_comb_name(grup, name_list)

        pred_lis = self.predict_calculation(data_grup[1:])
        pred_lis = np.concatenate(tuple(pred_lis), axis=1)
        self.adjust_combinationvalues(data_grup)

        for i in range(pred_lis.shape[1]):
            bi_list.append(np.expand_dims(np.subtract(data_grup[0].data_grup[:, i], pred_lis[:, i]), axis=1))
        bi_list = np.concatenate(tuple(bi_list), axis=1)

        return pred_lis, bi_list, name_list

    def get_comb_name(self, grup, name_list):
        counter = None

        for i in range(grup[0].concentration.shape[0]):
            fr = ''
            fr, counter = self.name_rec(fr, grup[1:], counter)

            name_list.append(fr)

        return name_list

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

            fr = str(group.name) + " " + str(group.concentration[counter[ind], 0]) + ' ' + str(group.unidade) + ' ' + fr

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

    def predict_calculation(self, data_grup):

        if len(data_grup) > 1:
            d1 = data_grup.pop(0).data_grup

            predicted = []

            d2 = self.predict_calculation(data_grup)

            i = 0
            while i < d1.shape[1]:

                for d in d2:
                    predicted.append(np.add(np.expand_dims(d1[:, i], axis=1),
                                            np.multiply(np.expand_dims(1 - d1[:, i], axis=1), d)))

                i += 1

            return predicted

        else:
            return [data_grup[0].data_grup]

    def selection(self, pred_lis, grup, data_grup, counter, ind=0):

        if len(grup) > 0:

            if counter is None:
                counter = np.zeros(len(data_grup), dtype=np.int16)

            group = grup.pop(0)

            if counter[ind] == group.concentration.shape[0]:
                counter[ind] = 0
                counter[ind - 1] += 1
            pred_lis, counter = self.selection(pred_lis, grup, data_grup, counter, ind + 1)

            if ind == len(counter) - 1:
                counter[ind] += 1

            return pred_lis, counter

        else:
            combp = counter
            pred_cal = self.bi_recursive(data_grup, combp, counter)

            pred_lis.append(np.expand_dims(pred_cal, axis=1))  # pred_cal)

            return pred_lis, counter

    def bi_recursive(self, data, combp, counter, ind=0):
        predict_calculation = ''
        if len(counter) > 2:

            count = np.delete(counter, 0)

            self.bi_recursive(data, combp, count, ind + 1)

        else:

            predict_calculation = np.add(data[ind].data_grup[:, counter[0]],
                                         np.multiply(data[ind + 1].data_grup[:, counter[1]],
                                                     (1 - data[ind].data_grup[:, counter[0]])))

        return predict_calculation

    @staticmethod
    def adjust_combinationvalues(data_grup):
        for i in range(len(data_grup[0].data_grup)):
            for k in range(len(data_grup[0].data_grup[i])):
                if data_grup[0].data_grup[i, k] < 0:
                    data_grup[0].data_grup[i, k] = 0
