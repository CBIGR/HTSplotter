import numpy as np
import sys
import os

class Readfile:
    def __init__(self, file_name):

        self.file_path = file_name
        # output
        self.information_file = None
        self.date = None
        self.elapsed = None
        self.header = None
        self.data = None
        self.data_std = None
        self.std = None
        self.std_type = None
        self.std_header = None
        self.date_info = None

        self.f = self.open_file()

        self.read_file()

        self.close_file()

        # self.get_split_data_std()

    def open_file(self):
        f = open(self.file_path, "r")
        return f

    def read_file(self):
        # get main information for each atribute
        n = 0
        informa_aux = []
        for line in iter(self.f.readline, ''):
            line_aux = line[:-1].split('\t')  # if is a txt file
            n += 1
            if line_aux[0] == 'Date Time':
                header_comple_row = line_aux
                for k in range(len(header_comple_row)):
                    header_comple_row[k] = header_comple_row[k].replace('"', '')
                break
                # indicate the header row number
            else:
                informa_aux.append(line_aux[0])

        initpos = self.f.tell()
        self.f.seek(initpos)

        rows = 0  # number of data rows

        for _ in iter(self.f.readline, ''):
            rows += 1

        self.information_file = informa_aux

        data_size = (len(header_comple_row) - 2)
        dt = []

        elap = np.zeros(rows, dtype=np.float32)
        data = np.zeros((rows, data_size), dtype=np.float32)

        self.header = header_comple_row[2:]

        self.f.seek(initpos)
        j = 0
        for line in iter(self.f.readline, ''):
            line_aux = line[:-1].split('\t')
            dt.append(line_aux[0])

            elap[j] = int(float(line_aux[1]))

            try:
                data[j, :] = np.asarray(list(map(float, line_aux[2:])))
            except ValueError:
                for sear in line_aux[2:]:
                    if len(sear) == 0:
                        sear = 0

            j += 1

        self.elapsed = elap
        self.date = dt
        self.date_info = [dt[0], dt[-1]]
        self.data_std = data


    def close_file(self):
        self.f.close()
        return

    # to be called after output file and error files
    def get_split_data_std(self, split):
        num_rows, num_cols = self.data_std.shape
        # split = int(num_cols/2)
        self.data = self.data_std[:, :split]
        self.std = self.data_std[:, split:]

    def get_data_nostd(self):
        self.data = self.data_std
        self.std = []

