import random
import numpy as np
from matplotlib.pyplot import cm

class ColorsHtsplots():
    def __init__(self, grup):

        n = int(np.ceil(grup**(1/3)))
        n1 = n
        n2 = n
        n3 = n
        while n1*n2*n3 > grup:
            if n3 < n1:
                n1 -= 1
            else:
                n3 -= 1
        if n1*n2*n3 < grup:
            n3 += 1
        test = list(cm.Set2(np.linspace(0, 1, 8)))
        if len(test) < grup:
            for i in range(n3):
                for j in range(n2):
                    for k in range(n1):
                        r = k / n3
                        b = i / n1
                        g = j / n2
                        test.append((r, b, g))
            test = test[1:-1]
            test = [test[-(1+i)] for i in range(len(test))]
        # file = open("C:/Users/cdcarval/OneDrive - UGent/Bureaublad/colorList.txt", "w")
        # file.write(str(test))
        # file.close()
        self.map = test