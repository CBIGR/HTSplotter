import numpy as np
import os
import math
from save_hdf5file import Hdf5database
from synergyprocessing import Commonprocessing
from scipy.optimize import fsolve, least_squares, curve_fit
import matplotlib.pyplot as plt
import itertools


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
        processing = Commonprocessing(grup)
        # correct inhibition below control

        # processing.adjust_data_new(data_grup[1:])

        name_list = processing.get_comb_name()

        pred_lis = self.predict_calculation(data_grup[1:])
        pred_lis = np.concatenate(tuple(pred_lis), axis=1)
        # processing.adjust_combinationvalues(data_grup)

        for i in range(pred_lis.shape[1]):
            bi_list.append(np.expand_dims(np.subtract(data_grup[0].data_grup[:, i], pred_lis[:, i]), axis=1))
        bi_list = np.concatenate(tuple(bi_list), axis=1)

        return pred_lis, bi_list, name_list

    # def get_comb_name(self, grup, name_list):
    #     counter = None
    #
    #     for i in range(grup[0].concentration.shape[0]):
    #         fr = ''
    #         fr, counter = self.name_rec(fr, grup[1:], counter)
    #
    #         name_list.append(fr)
    #
    #     return name_list
    #
    # def name_rec(self, fr, grup, counter, ind=0):
    #
    #     if len(grup) > ind:
    #         if counter is None:
    #             counter = np.zeros(len(grup), dtype=np.int16)
    #
    #         group = grup[ind]
    #
    #         for i in range(ind, -1, -1):
    #             if counter[i] == grup[i].concentration.shape[0]:
    #                 counter[i] = 0
    #                 if i > 0:
    #                     counter[i - 1] += 1
    #
    #         fr, counter = self.name_rec(fr, grup, counter, ind + 1)
    #
    #         fr = str(group.name) + " " + str(group.concentration[counter[ind], 0]) + ' ' + str(group.unidade) + ' ' + fr
    #
    #         if ind == len(counter) - 1:
    #             counter[ind] += 1
    #
    #         return fr, counter
    #
    #     else:
    #
    #         return fr, counter

    # @staticmethod
    # def adjust_data_new(data):
    #
    #     for i in range(len(data)):
    #         for j in range(len(data[i].data_grup)):
    #             for k in range(len(data[i].data_grup[j])):
    #                 if data[i].data_grup[j, k] < 0:
    #                     data[i].data_grup[j, k] = 0
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

    # def selection(self, pred_lis, grup, data_grup, counter, ind=0):
    #
    #     if len(grup) > 0:
    #
    #         if counter is None:
    #             counter = np.zeros(len(data_grup), dtype=np.int16)
    #
    #         group = grup.pop(0)
    #
    #         if counter[ind] == group.concentration.shape[0]:
    #             counter[ind] = 0
    #             counter[ind - 1] += 1
    #         pred_lis, counter = self.selection(pred_lis, grup, data_grup, counter, ind + 1)
    #
    #         if ind == len(counter) - 1:
    #             counter[ind] += 1
    #
    #         return pred_lis, counter
    #
    #     else:
    #         combp = counter
    #         pred_cal = self.bi_recursive(data_grup, combp, counter)
    #
    #         pred_lis.append(np.expand_dims(pred_cal, axis=1))  # pred_cal)
    #
    #         return pred_lis, counter
    #
    # def bi_recursive(self, data, combp, counter, ind=0):
    #     predict_calculation = ''
    #     if len(counter) > 2:
    #
    #         count = np.delete(counter, 0)
    #
    #         self.bi_recursive(data, combp, count, ind + 1)
    #
    #     else:
    #
    #         predict_calculation = np.add(data[ind].data_grup[:, counter[0]],
    #                                      np.multiply(data[ind + 1].data_grup[:, counter[1]],
    #                                                  (1 - data[ind].data_grup[:, counter[0]])))
    #
    #     return predict_calculation

    # @staticmethod
    # def adjust_combinationvalues(data_grup):
    #     for i in range(len(data_grup[0].data_grup)):
    #         for k in range(len(data_grup[0].data_grup[i])):
    #             if data_grup[0].data_grup[i, k] < 0:
    #                 data_grup[0].data_grup[i, k] = 0


class Hsamethod:
    def __init__(self, grup, data_grup, i, j):
        self.grup = grup
        self.grup_dataalone = data_grup[1:]
        self.data_grup = data_grup
        self.i = i
        self.j = j

        self.hsavalues = []
        self.maximumeffect = []
        self.combname = []

        self.hsa_score()

    def hsa_score(self):

        maximumeffect_lis, hsa_list, name_list = self.hsa_score_matrix(self.grup, self.data_grup)

        self.hsavalues=np.asarray(hsa_list)
        self.maximumeffect = np.asarray(maximumeffect_lis)
        self.combname = name_list

    def hsa_score_matrix(self, grup, data_grup):
        hsa_list = []
        processing = Commonprocessing(grup)
        # correct inhibition below control

        processing.adjust_data_new(self.grup_dataalone)

        name_list = processing.get_comb_name()

        maximumeffect_lis = self.get_maximumeffect(data_grup[1:])
        maximumeffect_lis = np.concatenate(tuple(maximumeffect_lis), axis=1)
        processing.adjust_combinationvalues(data_grup)

        for i in range(maximumeffect_lis.shape[1]):
            hsa_list.append(np.expand_dims(np.subtract(data_grup[0].data_grup[:, i], maximumeffect_lis[:, i]), axis=1))
        hsa_list = np.concatenate(tuple(hsa_list), axis=1)

        return maximumeffect_lis, hsa_list, name_list

    def get_maximumeffect(self, data_grup):
        # first get the maximum
        if len(data_grup) > 1:
            d1 = data_grup.pop(0).data_grup

            maximum = []

            d2 = self.get_maximumeffect(data_grup)

            i = 0
            while i < d1.shape[1]:

                for d in d2:
                    maximum.append(np.maximum(np.expand_dims(d1[:, i], axis=1), d))

                i += 1

            return maximum

        else:
            return [data_grup[0].data_grup]


class LoeweMethod:

    def __init__(self, combination, grup, data_grup, i, j):
        self.grup = grup
        self.data_grup = data_grup
        self.time = combination.elapse
        self.curvescompoundalone = combination.curvessynergy
        self.curvescompoundalonename = combination.doseresponsecurve.compoundname
        self.concentration = combination.curvessynergconcentration
        self.time_selected = combination.time_selected
        self.time_position = combination.time_position  # combination.time_position
        self.curvetime = [0, 1, 2]

        self.yLoewe = None
        self.yc = None
        self.CI = None
        self.CIscore = None
        self.yloewescore = []

        processing = Commonprocessing(self.grup)
        # correct inhibition below control

        processing.adjust_data_new(self.data_grup[1:])

        # name_list = processing.get_comb_name()
        name_list = []
        yloewe = []
        difference = []
        CI = []

        # for time in range(len(self.time)):
            # for each compound alone, each concentration I need to find yLoewe
        aux = self.grup[1:]
        indexes = [0 for _ in range(len(aux))]
        max_indexes = [drug.concentration.shape[0] for drug in aux]
        a = 0
        tot = np.product(max_indexes)
        # adjust the values
        # get the Yobsered
        while a < tot:
            yloewe.append([])
            CI.append([])
            difference.append([])
            xs = [aux[k].concentration[indexes[k], 0] for k in range(len(aux))]
            name_list.append('_'.join([aux[k].name + ' ' + str(aux[k].concentration[indexes[k], 0]) for k in range(len(aux))]))
            for time in range(len(self.time_selected)):
                # for combs in range(len(self.grup[0].concentration)):
                #     if np.array_equal(self.grup[0].concentration[combs][:-1], np.asarray(xs)):
                #         yc = self.data_grup[0][int(self.grup[0].concentration[a][1])]
                # yc = self.data_grup[0].data_grup[time][a]
                yloewe[a].append(self.get_yLoewe(xs, self.curvescompoundalone, self.curvetime[time]))
                CI[a].append(self.get_CI(self.data_grup[0].data_grup[self.time_position[time]][a], xs,
                                         self.curvescompoundalone, self.time_position[time]))
                difference[a].append(self.data_grup[0].data_grup[self.time_position[time]][a] -
                                     self.get_yLoewe(xs, self.curvescompoundalone, self.time_position[time]))
                # print(a, self.time_position[time])
                # loewescore[a].append(self.data_grup[0].data_grup[self.time_position[time]][a]
                #                                                  - yloewe[-1][0])
            indexes[-1] += 1
            for s in range(len(indexes)-1, -1, -1):
                if indexes[s] >= max_indexes[s]:
                    indexes[s] = 0
                    indexes[s-1] += 1
            a += 1
        # for eachloewe in range(len(yloewe)):
        #     loewescore.append([])

            # for tim in range(len(self.time_selected)):
                # print(eachloewe, self.time_position[tim])
                # loewescore[eachloewe].append(self.data_grup[0].data_grup[self.time_position[tim]][eachloewe]
                #                  - yloewe[eachloewe][tim])
                # print(loewescore)

        self.yloewe = np.asarray(yloewe).T
        self.CIscore = np.asarray(CI).T
        self.difference = np.asarray(difference).T
        self.combname = name_list
        # print('we')

    def get_term(self, x, xi, curve, case=False):

        # xi = ith drug concentration

        if case:

            func = xi/(curve.m*(((x-curve.Emin)/(curve.Emax - x))**(1/curve.lambd)))

        else:
            func = xi/(curve.m*(((x-curve.Emin)/(curve.Emax - self.yLoewe))**(1/curve.lambd)))

        return func

    def get_yLoewe(self, xi, curves, time):

        dy = 0.01
        print(curves[0][time].time_point)
        Emax = np.max([curve[time].Emax for curve in curves if curve[time].Emax > 0.0])
        Emin = np.min([curve[time].Emin for curve in curves if curve[time].Emin > 0.0])
        print('from the same curve??', Emax, Emin)
        Vmax = Emax
        Vmin = Emin
        print('from the same curve??', Vmax, Vmin)
        fy = 1e5

        while abs(fy - 1) > 1e-7:

            y = np.linspace(Vmax, Vmin, 100)
            fy = np.asarray([self.get_solver_terms(x, xi, curves, time)-1 for x in y])

            i = 0
            while fy[i]*fy[i+1] > 0 and i < len(fy)-1:
                i += 1

            Vmax = max([y[i], y[i + 1]])
            Vmin = min([y[i], y[i + 1]])

            y = (Vmax+Vmin)/2
            print('before', fy, y, xi)
            fy = self.get_solver_terms(y, xi, curves, time)
            print('after', fy, y, xi)
        print(curves[0][time].time_point, fy)

        # validate yloewe has to be equal 1

        self.yLoewe = y

        return y


    def get_solver_terms(self, y, xi, curves, time):

        res = []
        a = 0
        for curve in curves:
            # print(curve[time].compound_alone)
            ct = curve[time]
            res.append(self.get_term(y, xi[a], ct, True))
            a += 1
        res = [0 if math.isnan(x) else x for x in res]
        return sum(res)


    def get_CI(self,y , xi, curves, time):

        res = []
        a = 0
        for curve in curves:
            ct = curve[time]
            res.append(self.get_term(y, xi[a], ct))
            a += 1

        return sum(res)


class ZipMethomultidimen:

    def __init__(self, combination, grup, data_grup):

        self.grup = grup
        self.data_grup = data_grup
        self.time_selected = combination.elapse
        self.emin0 = 0
        self.zipvalues = []
        self.yzipvalues = []
        self.combname = []

        for time in range(len(self.time_selected)):

            cols = {drug: col for drug, col, in zip(self.grup[0].name.split('_'),
                                                    range(len(self.grup[0].name.split('_'))))}

            for i in range(len(self.grup)):

                drugs = self.grup[i].name.split('_')
                aux = self.grup[i].concentration[:, :-1]
                auxy = self.data_grup[i].data_grup[time].flatten()

                if aux.shape[1] < len(cols):

                    aux1 = np.zeros((aux.shape[0], len(cols)))

                    k = 0
                    for drug in drugs:
                        aux1[:, cols[drug]] = aux[:, k]
                        k += 1

                    concentrations = np.r_[concentrations, aux1]
                    response = np.r_[response, auxy]

                else:
                    concentrations = aux.copy()
                    response = auxy.copy()

            total = np.c_[concentrations, response]
            del response, concentrations
            ft = total.shape[1] * 'float, '
            ft = ft[:-2]
            order = ['f' + str(k) for k in range(total.shape[1])]
            total.view(ft).sort(kind='quicksort', axis=0, order=order)

            concentrations = []
            aux = np.asarray([0])
            for i in range(len(cols)):
                concentrations.append(np.r_[aux, self.grup[-1*(i+1)].concentration[:, :-1].flatten()])

            indexes = itertools.product(*[list(range(concentration.size)) for concentration in concentrations])
            response = np.zeros(tuple([concentration.size for concentration in concentrations])[::-1])
            aux_concentrations = itertools.product(*concentrations)
            concentrations_array = [np.zeros(response.shape) for _ in range(len(cols))]
            for index, concentration in zip(indexes, aux_concentrations):
                if np.sum(concentration) > 0:
                    mask = ~np.sum(total[:, -2::-1] - concentration, axis=1, dtype=bool)
                    response[index[::-1]] = total[mask, -1]
                    for i in range(len(cols)):
                        concentrations_array[i][index[::-1]] = concentration[-1*(i+1)]

            self.m = {}
            self.lambd = {}
            self.emax = {}
            self.emin = {}

            # fit para as drogas sozinhas
            emin0 = None
            for cycle in range(2):

                m0 = []
                lamb0 = []
                emax0 = []

                if cycle == 1:
                    emin0 = np.mean(emin0)
                else:
                    emin0 = []

                for drug in cols.keys():

                    slc = [0 for _ in range(len(cols))]
                    slc[cols[drug]] = slice(1, None)

                    drug_concentration = (concentrations_array[cols[drug]][tuple(slc)]).flatten()
                    drug_response = (response[tuple(slc)]).flatten()

                    if cycle == 0:
                        popt, _ = curve_fit(self.curve_alone_first,
                                            drug_concentration, drug_response, method='lm', maxfev=int(1e7))
                    else:
                        popt, _ = curve_fit(self.curve_alone,
                                            drug_concentration, drug_response, method='trf', maxfev=int(1e7),
                                            bounds=((0, 0, 0), (np.inf, np.inf, np.inf)))


                    if cycle == 0:
                        m0.append(popt[0])
                        lamb0.append(popt[1])
                        emax0.append(popt[2])
                        emin0.append(popt[3])
                    else:
                        self.m[str(cols[drug])] = np.asarray([popt[0]])
                        self.lambd[str(cols[drug])] = np.asarray([popt[1]])
                        self.emax[str(cols[drug])] = np.asarray([popt[2]])
                        self.emin[str(cols[drug])] = np.asarray([emin0])

            indexes = np.asarray(list(itertools.product(*[list(range(dim)) for dim in response.shape])))
            mask = np.zeros(indexes.shape[0], dtype=bool)
            for i in range(indexes.shape[0]):
                if 0 < np.sum(indexes[i, :] == 0) < len(response.shape):
                    mask[i] = True
            indexes = indexes[mask]

            slices = []
            for i in range(indexes.shape[0]):
                idx = [j == 0 for j in indexes[i, :]]
                idx = [j for j in range(indexes[i, :].size) if idx[j]]
                for k in idx:
                    aux = list(indexes[i, :])
                    aux[k] = slice(1, None)
                    slices.append(aux)

            # fit para as drogas em combinacoes
            concentrations = concentrations[::-1]
            for i in range(len(slices[0])-2, -1, -1):
                indexes = [np.sum([s == 0 for s in slc]) == i for slc in slices]
                for j in range(len(indexes)):
                    if indexes[j]:
                        k = 0
                        idx = None
                        while k < len(slices[j]) and idx is None:
                            if type(slices[j][k]) is slice:
                                idx = k
                            k += 1
                        drug_concentration = concentrations[idx][1:]
                        drug_response = response[tuple(slices[j])]

                        idx_ = []
                        concentrations_ = []
                        for k in range(len(slices[j])):
                            if type(slices[j][k]) is slice:
                                idx_.append(-1)
                                concentrations_.append(-1)
                            else:
                                idx_.append(slices[j][k])
                                concentrations_.append(concentrations[k][slices[j][k]])
                        self.curve_combine_recursive(np.asarray(idx_), slices[j], np.asarray(concentrations_),
                                                     drug_concentration, drug_response,
                                                     [len(k) for k in concentrations], idx)

            indexes = list(itertools.product(*[[l for l in range(a)] for a in [len(k) for k in concentrations]]))
            concentrations_ = [[concentrations[l][i] for l, i in zip(range(len(idx)), idx)] for idx in indexes]
            self.y_zip = np.zeros(response.shape)
            self.get_zip(indexes, concentrations_)
            self.y_prediction = np.zeros(response.shape)
            self.get_prediction(indexes, concentrations_)

            self.delta = self.y_prediction - self.y_zip

            self.delta_list = []
            self.y_predictionlist = []
            namelist = []
            indexes = list(itertools.product(*[[k for k in range(len(concentrations[i]))]
                                               for i in range(len(cols))]))
            name = [values for values in cols.items()]
            for index in indexes:
                aux = [concentrations[k][index[k]] for k in range(len(concentrations))]
                # namelist.append(aux)
                if 0.0 not in aux:
                    nameconcentration = [name[na][0] + ' ' + str(aux[na]) for na in range(len(aux))]
                    namelist.append('_'.join(nameconcentration))
                    if self.delta[index] >= 1:
                        self.delta_list.append(0.0)
                    else:
                        self.delta_list.append(self.delta[index])
                    self.y_predictionlist.append(self.y_prediction[index])
                    print('no', time)
                # aux.append(self.delta[index])

            self.zipvalues.append(self.delta_list)
            self.yzipvalues.append(self.y_predictionlist)
            self.combname = namelist
        print('done')

    def get_combinationscurves(self):
        print('ojnono')


    def get_Emin(self, leng, time):
        emin = []
        position = 0
        while position < leng:
            emin.append(self.curvessynergynormal[position][time].Emin)
            # print('sto', emin)
            position +=1

        return sum(emin)/leng

    # def curve_alone(self, x, m, l, emax):
    #
    #     a = (x/m)**l
    #
    #     return (self.Eminaverage[-1] + emax*a)/(1+a)

    def curve_alone_first(self, x, m, l, emax, emin):
        a = (x / m) ** l

        return (emin + emax * a) / (1 + a)

    def curve_alone(self, x, m, l, emax):

        a = (x/m)**l

        return (self.emin0 + emax*a)/(1+a)

    def curve_combine(self, x, m, l, emax):

        a = (x / m) ** l

        return (self.y_ + emax*a)/(1+a)

    def curve_combine_recursive(self, idx, concentrations_idx, concentrations, x, y, dims=None, sidx=None, y_=None,
                                order=''):

        if idx.size == 0:

            self.y_ = y_
            popt, _ = curve_fit(self.curve_combine,
                                x, y, method='trf', maxfev=int(1e7),
                                bounds=((0, 0, 0), (np.inf, np.inf, np.inf)))

            i = [concentrations_idx[int(k)] for k in order[1:][::-1]]
            self.m[order][tuple(i)] = popt[0]
            self.lambd[order][tuple(i)] = popt[1]
            self.emax[order][tuple(i)] = popt[2]
            self.emin[order][tuple(i)] = y_

        else:

            aux0 = [i for i in range(len(idx)) if idx[i] > 0]

            combs = itertools.permutations(aux0, len(aux0))
            keys = list(self.m.keys())
            for comb in combs:
                prev_drug_comb = ''
                for a in comb:
                    prev_drug_comb = str(a) + prev_drug_comb
                drug_comb = str(sidx) + prev_drug_comb
                if drug_comb not in keys:
                    shape = tuple([dims[int(k)] for k in prev_drug_comb[::-1]])
                    self.m[drug_comb] = np.zeros(shape)
                    self.lambd[drug_comb] = np.zeros(shape)
                    self.emax[drug_comb] = np.zeros(shape)
                    self.emin[drug_comb] = np.zeros(shape)

                if len(prev_drug_comb) > 1:
                    drugs = [concentrations_idx[int(k)] for k in prev_drug_comb[1:]]
                else:
                    drugs = [0]

                concentration = concentrations[int(prev_drug_comb[0])]

                y_ = self.curve_alone_first(concentration, self.m[prev_drug_comb][tuple(drugs)],
                                            self.lambd[prev_drug_comb][tuple(drugs)],
                                            self.emax[prev_drug_comb][tuple(drugs)],
                                            self.emin[prev_drug_comb][tuple(drugs)])

                self.curve_combine_recursive(np.asarray([]), concentrations_idx, None, x, y, y_=y_, order=drug_comb)

    def get_zip(self, indexes, concentrations):

        m = [self.m[str(a)][0] for a in range(len(indexes[0]))]
        lambd = [self.lambd[str(a)][0] for a in range(len(indexes[0]))]
        emax = [self.emax[str(a)][0] for a in range(len(indexes[0]))]
        emin = [self.emin[str(a)][0] for a in range(len(indexes[0]))]

        idx = ''
        for i in range(len(indexes[0])):
            idx += str(i)

        for i in range(len(indexes)):
            for j in range(len(indexes[i])):
                comb = itertools.combinations(idx, j+1)
                for c in comb:
                    aux = 1
                    for k in c:
                        ik = int(k)
                        a = (concentrations[i][ik]/m[ik])**lambd[ik]
                        aux *= (emin[ik] + emax[ik]*a)/(1+a)

                    if len(c) % 2 == 0:
                        self.y_zip[tuple(indexes[i])] -= aux
                    else:
                        self.y_zip[tuple(indexes[i])] += aux

    def get_prediction(self, indexes, concentrations):

        idx = ''
        for i in range(len(indexes[0])):
            idx += str(i)

        for i in range(len(indexes)):
            comb = list(itertools.permutations(idx, len(idx)))
            aux = 0
            for c in comb:

                nc = ''
                for b in c:
                    nc += b

                iaux = tuple([indexes[i][int(k)] for k in nc[1:][::-1]])
                m = self.m[nc][iaux]
                lambd = self.lambd[nc][iaux]
                emax = self.emax[nc][iaux]
                emin = self.emin[nc][iaux]

                a = (concentrations[i][int(c[0])]/m)**lambd
                aux += (emin + emax*a)/(1+a)

            self.y_prediction[tuple(indexes[i])] += aux/len(comb)


    def ZIPref(self, x11, x22, m1, m2, l1, l2, emax1, emax2, emin1, emin2):

        a_ = (x11/m1)**l1
        b_ = (x22/m2)**l2

        a = (emin1 + emax1*a_)/(1 + a_)
        b = (emin2 + emax2*b_)/(1 + b_)

        return a + b - a*b

    def ZIPFit(self, x11, x22, m1, m2, l1, l2, emax1, emax2, emin1, emin2, m12, m21, l12, l21, emax12, emax21):

        f = np.zeros(x11.shape)
        g = np.zeros(x11.shape)

        for i in range(x11.shape[0]):
            for j in range(x11.shape[1]):

                a = (x11[i, j]/m1)**l1
                b = (x22[i, j]/m2)**l2

                y1 = (emin1 + emax1*a)/(1 + a)
                y2 = (emin2 + emax2*b)/(1 + b)

                c = (x11[i, j]/m12[i])**l12[i]
                d = (x22[i, j]/m21[j])**l21[j]

                f[i, j] = (y2 + emax12[i]*c)/(1+c)
                g[i, j] = (y1 + emax21[j]*d)/(1+d)

        return 0.5*(g + f)



##########################################
#### old ZIP
class ZipMethod:

    def __init__(self, combination, grup, data_grup):
        self.grup = grup
        self.data_grup = data_grup
        # self.time = combination.elapse
        self.curvessynergynormal = combination.curvessynergynormal
        self.curvescompoundalonename = combination.doseresponsecurve.compoundname
        self.compoundleng = np.copy(combination.curvessynergynormal)
        self.concentration = combination.curvessynergconcentration
        self.time_selected = combination.elapse

        self.yfromzipfit = []
        self.yzipreference = []
        self.zipscore = []
        self.combinationname = None  # to check on heatmap the combination results

        self.zipfit = []
        self.zipfit = []
        self.Eminaverage = []
        self.Emineach = []

        namelist = []
        x1_name = self.grup[1].name
        x2_name = self.grup[2].name
        x1 = np.asarray(list(self.concentration[0]))
        x2 = np.asarray(list(self.concentration[1]))
        x11, x22 = np.meshgrid(x1, x2)

        # aux_list = []
        # from the dose response curve, obtained we need to compute Emin avarage
        # and then a new curve for monoterapy
        for time in range(len(self.time_selected)):
            x11, x22 = np.meshgrid(x1, x2)
            print(time)

            m1 = self.curvessynergynormal[0][time].m
            m2 = self.curvessynergynormal[1][time].m
            l1 = self.curvessynergynormal[0][time].lambd
            l2 = self.curvessynergynormal[1][time].lambd
            emax1 = self.curvessynergynormal[0][time].Emax
            emax2 = self.curvessynergynormal[1][time].Emax
            ###

            y = np.zeros(x11.shape)

            y[1:, 1:] = self.data_grup[0].data_grup[time, :].reshape(y[1:, 1:].shape).T
            y[0, 1:] = self.data_grup[1].data_grup[time, :]
            y[1:, 0] = self.data_grup[2].data_grup[time, :]
            self.Eminaverage = (self.curvessynergynormal[0][time].Emin + self.curvessynergynormal[1][time].Emin) / 2

            # compute a new dose response curve, based on the new minin
            popt, _ = curve_fit(self.curve_alone, x1[1:], y[0, 1:], method='lm', maxfev=int(1e7))
            m1, l1, emax1 = popt
            popt, _ = curve_fit(self.curve_alone, x2[1:], y[1:, 0], method='lm', maxfev=int(1e7))
            m2, l2, emax2 = popt
            emin1 = self.Eminaverage
            emin2 = self.Eminaverage

            # recursive
            m12 = [m1]
            m21 = [m2]
            l12 = [l1]
            l21 = [l2]
            emax12 = [emax1]
            emax21 = [emax2]

            j = 1
            for x1_ in x1[1:]:
                self.y_ = (emin1 + emax1 * (x1_ / m1) ** l1) / (1 + (x1_ / m1) ** l1)
                popt, _ = curve_fit(self.curve_combine, x2[1:], y[1:, j], method='lm', maxfev=int(1e7))
                # ,maxfev=10000 method='trf',
                # bounds=((np.min(x2[y[:, j] > 0.1*self.y_]), 0.0, np.mean(y[:, j])),
                #         (np.inf, np.inf, 1)), max_nfev=100000)
                m21.append(popt[0])
                l21.append(popt[1])
                emax21.append(popt[2])
                j += 1

            i = 1
            for x2_ in x2[1:]:
                self.y_ = (emin2 + emax2 * (x2_ / m2) ** l2) / (1 + (x2_ / m2) ** l2)
                popt, _ = curve_fit(self.curve_combine, x1[1:], y[i, 1:], method='lm',
                                    maxfev=int(1e7))  # , method='trf',
                # bounds=((np.min(x1[y[i, :] > 0.1*self.y_]), 0.0, np.mean(y[i, :])),
                #         (np.inf, np.inf, 1)), max_nfev=100000)
                m12.append(popt[0])
                l12.append(popt[1])
                emax12.append(popt[2])
                i += 1

            zipref = self.ZIPref(x11, x22, m1, m2, l1, l2, emax1, emax2, emin1, emin2)
            zipfit = self.ZIPFit(x11, x22, m1, m2, l1, l2, emax1, emax2, emin1, emin2,
                                 m12, m21, l12, l21, emax12, emax21)

            auxlist1 = zipfit - zipref

            x11 = x11[1:, 1:].flatten()
            x22 = x22[1:, 1:].flatten()

            namelist.append([x1_name + ' ' + str(x11[i]) + '_' + x2_name + ' ' + str(x22[i]) for i in range(len(x11))])
            self.zipscore.append(auxlist1[1:, 1:].flatten())
            self.yfromzipfit.append(zipfit[1:, 1:].flatten())
            self.yzipreference.append(zipref[1:, 1:].flatten())
            # self.zipref.append(zipref[1:, 1:].flatten())
            # self.zipfit.append(zipfit[1:, 1:].flatten())
            pass

        self.curvefromnewemin = []  # new fit for each compound alone, based on the new Emin.
        self.observationmatrixwithcompoundalone = []

        # at the end, we need to get yfromzipfit, yzipreference and zipscore
        self.yfromzipfit.append(self.zipfit[1:, 1:].flatten())
        self.yzipreference.append(self.zipref[1:, 1:].flatten())
        self.zipscore.append(aux_list)
        self.combinationname = namelist[0]  # to check on heatmap the combination results

    def get_combinationscurves(self):
        print('ojnono')

    def get_Emin(self, leng, time):
        emin = []
        position = 0
        while position < leng:
            emin.append(self.curvessynergynormal[position][time].Emin)
            # print('sto', emin)
            position += 1

        return sum(emin) / leng

    def curve_alone(self, x, m, l, emax):

        a = (x / m) ** l

        return (self.Eminaverage + emax * a) / (1 + a)

    def curve_combine(self, x, m, l, emax):

        a = (x / m) ** l

        return (self.y_ + emax * a) / (1 + a)

    def ZIPref(self, x11, x22, m1, m2, l1, l2, emax1, emax2, emin1, emin2):

        a_ = (x11 / m1) ** l1
        b_ = (x22 / m2) ** l2

        a = (emin1 + emax1 * a_) / (1 + a_)
        b = (emin2 + emax2 * b_) / (1 + b_)

        return a + b - a * b

    def ZIPFit(self, x11, x22, m1, m2, l1, l2, emax1, emax2, emin1, emin2, m12, m21, l12, l21, emax12, emax21):

        f = np.zeros(x11.shape)
        g = np.zeros(x11.shape)

        for i in range(x11.shape[0]):
            for j in range(x11.shape[1]):
                a = (x11[i, j] / m1) ** l1
                b = (x22[i, j] / m2) ** l2

                y1 = (emin1 + emax1 * a) / (1 + a)
                y2 = (emin2 + emax2 * b) / (1 + b)

                c = (x11[i, j] / m12[i]) ** l12[i]
                d = (x22[i, j] / m21[j]) ** l21[j]

                f[i, j] = (y2 + emax12[i] * c) / (1 + c)
                g[i, j] = (y1 + emax21[j] * d) / (1 + d)

        return 0.5 * (g + f)


class CombinationResponseCurve:

    def __init__(self, l12, l21, m12, m21, time):
        self.lambd12 = l12
        self.lambd21 = l21
        self.m12 = m12
        self.m21 = m21

        self.time_point = time
