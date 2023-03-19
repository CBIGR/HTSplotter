import numpy as np
np.set_printoptions(suppress=True)


class Groupping:

    def __init__(self, cab, compound, cellin, seeding, condition):

        self.name = cab[0][3]
        self.cell_name = None
        self.seeding = cab[0][1]
        self.concentration = None

        aux = []
        for c in cab:

            if cellin == c[0]:
                self.cell_name = c[0]
                a = c[4].split(" ")
                if "_" in a[0]:
                    b = a[0].split("_")
                    self.unidade = a[-1]

                    aux.append([float(b[i]) for i in range(len(b))])
                    aux[-1].append(c[-1])
                else:
                    self.unidade = a[-1]
                    aux.append([float(a[2 * i]) for i in range(int(len(a) / 2))])
                    aux[-1].append(c[-1])

        if len(aux) > 0:
            aux = np.asarray(aux)

            fmt = ""
            orde = ["f" + str(i) for i in range(aux.shape[1] - 1)]

            for i in range(aux.shape[1]):
                fmt += "float, "

            if aux.shape[1] > 1:
                aux.view(fmt[:-2]).sort(order=orde, axis=0)

            self.concentration = aux

class Grouppingmedium:

    def __init__(self, cab, compound, cellin, seeding, condition):

        self.name = cab[3]
        self.cell_name = None
        self.seeding = cab[1]
        self.concentration = None

        aux = []

        if cellin == cab[0]:
            self.cell_name = cab[0]
            a = cab[4].split(" ")
            if "_" in a[0]:
                b = a[0].split("_")
                self.unidade = a[-1]

                aux.append([float(b[i]) for i in range(int(len(b)))])
                aux[-1].append(cab[-1])
            else:
                self.unidade = a[-1]
                aux.append([float(a[2 * i]) for i in range(int(len(a) / 2))])
                aux[-1].append(cab[-1])

        aux = np.asarray(aux)

        fmt = ""
        orde = ["f" + str(i) for i in range(aux.shape[1] - 1)]

        for i in range(aux.shape[1]):
            fmt += "float, "

        if aux.shape[1] > 1:
            aux.view(fmt[:-2]).sort(order=orde, axis=0)

        self.concentration = aux

class Data_group():

    def __init__(self, data, position):
        self.data_grup = np.asarray([data[int(i[-1])] for i in position]).T


class Grouppingpertub:

    def __init__(self, cab, cellin, seeding, condition):

        self.name = cab[0][3]
        self.cell_name = None
        self.seeding = cab[0][1]
        self.concentration = None

        aux = []
        for c in cab:

            if cellin == c[0]:
                self.cell_name = c[0]
                a = c[4].split(" ")
                if "_" in a[0]:
                    b = a[0].split("_")
                    self.unidade = a[-1]

                    aux.append([float(b[i]) for i in range(len(b))])
                    aux[-1].append(c[-1])
                else:
                    self.unidade = a[-1]
                    aux.append([float(a[2 * i]) for i in range(int(len(a) / 2))])
                    aux[-1].append(c[-1])

        if len(aux) > 0:
            aux = np.asarray(aux)

            fmt = ""
            orde = ["f" + str(i) for i in range(aux.shape[1] - 1)]

            for i in range(aux.shape[1]):
                fmt += "float, "

            if aux.shape[1] > 1:
                aux.view(fmt[:-2]).sort(order=orde, axis=0)

            self.concentration = aux


class Grouppingzip:
    def __init__(self, cab, header, cellin, seeding, condition):
        self.headergroupping = []
        self.name = cab
        self.cell_name = cellin
        self.condition = condition[0]
        self.seeding = seeding[0]
        self.concentration = None
        self.unidade = None
        aux = []
        for h in header:
            if h[0] == self.cell_name:
                if h[1] == self.seeding:
                    if h[2] == self.condition:
                        if h[3] == self.name:
                            self.headergroupping.append(h)
                            a = h[4].split(" ")
                            print(h)
                            if "_" in a[0]:
                                b = a[0].split("_")
                                self.unidade = a[-1]

                                aux.append([float(b[i]) for i in range(len(b))])
                                aux[-1].append(h[-1])
                            else:
                                self.unidade = a[-1]
                                aux.append([float(a[2 * i]) for i in range(int(len(a) / 2))])
                                aux[-1].append(h[-1])
                            # concen = i[-1].split('_')
        if len(aux) > 0:
            aux = np.asarray(aux)

            fmt = ""
            orde = ["f" + str(i) for i in range(aux.shape[1] - 1)]

            for i in range(aux.shape[1]):
                fmt += "float, "

            if aux.shape[1] > 1:
                aux.view(fmt[:-2]).sort(order=orde, axis=0)

            self.concentration = aux

        #
        # aux = []
        # for c in cab:
        #
        #     if cellin == c[0]:
        #         self.cell_name = c[0]
        #         a = c[4].split(" ")
        #         if "_" in a[0]:
        #             b = a[0].split("_")
        #             self.unidade = a[-1]
        #
        #             aux.append([float(b[i]) for i in range(len(b))])
        #             aux[-1].append(c[-1])
        #         else:
        #             self.unidade = a[-1]
        #             aux.append([float(a[2 * i]) for i in range(int(len(a) / 2))])
        #             aux[-1].append(c[-1])
        #
        # if len(aux) > 0:
        #     aux = np.asarray(aux)
        #
        #     fmt = ""
        #     orde = ["f" + str(i) for i in range(aux.shape[1] - 1)]
        #
        #     for i in range(aux.shape[1]):
        #         fmt += "float, "
        #
        #     if aux.shape[1] > 1:
        #         aux.view(fmt[:-2]).sort(order=orde, axis=0)
        #
        #     self.concentration = aux


class GrouppingGrowthRate:

    def __init__(self, conditiondrug, header, cellin, condition, seeding):
        self.headergroupping = []
        self.name = conditiondrug
        self.cell_name = cellin
        self.condition = condition
        self.seeding = seeding
        self.concentration = None
        self.unidade = None
        aux = []
        for h in header:
            if h[0] == self.cell_name:
                if h[1] == self.seeding:
                    if h[2] == self.condition:
                        if h[3] == self.name:
                            self.headergroupping.append(h)
                            a = h[4].split(" ")
                            print(h)
                            if "_" in a[0]:
                                b = a[0].split("_")
                                self.unidade = a[-1]

                                aux.append([float(b[i]) for i in range(len(b))])
                                aux[-1].append(h[-1])
                            else:
                                self.unidade = a[-1]
                                aux.append([float(a[2 * i]) for i in range(int(len(a) / 2))])
                                aux[-1].append(h[-1])
                            # concen = i[-1].split('_')
        if len(aux) > 0:
            aux = np.asarray(aux)

            fmt = ""
            orde = ["f" + str(i) for i in range(aux.shape[1] - 1)]

            for i in range(aux.shape[1]):
                fmt += "float, "

            if aux.shape[1] > 1:
                aux.view(fmt[:-2]).sort(order=orde, axis=0)

            self.concentration = aux

class GrouppingGrowthRateCombination:

    def __init__(self, conditiondrug, header, cellin, seeding):
        self.headergroupping = []
        self.name = conditiondrug
        self.cell_name = cellin
        self.condition = None
        self.seeding = seeding
        self.concentration = None
        self.unidade = None
        aux = []
        for h in header:
            if h[0] == self.cell_name:
                if h[1] == self.seeding:
                    # if h[2] == self.condition:
                    if h[3] == self.name:
                        self.headergroupping.append(h)
                        self.condition = h[2]
                        a = h[4].split(" ")
                        print(h)
                        if "_" in a[0]:
                            b = a[0].split("_")
                            self.unidade = a[-1]

                            aux.append([float(b[i]) for i in range(len(b))])
                            aux[-1].append(h[-1])
                        else:
                            self.unidade = a[-1]
                            aux.append([float(a[2 * i]) for i in range(int(len(a) / 2))])
                            aux[-1].append(h[-1])
                            # concen = i[-1].split('_')

        if len(aux) > 0:
            aux = np.asarray(aux)

            fmt = ""
            orde = ["f" + str(i) for i in range(aux.shape[1] - 1)]

            for i in range(aux.shape[1]):
                fmt += "float, "

            if aux.shape[1] > 1:
                aux.view(fmt[:-2]).sort(order=orde, axis=0)

            self.concentration = aux