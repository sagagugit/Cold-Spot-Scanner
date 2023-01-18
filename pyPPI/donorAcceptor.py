import os
import re

"""
Reads donor acceptor matrix
"""


class DonAcceptor:
    _instance = None

    def __init__(self, path=os.path.join(os.path.dirname(__file__), 'DonAcc2.txt')):
        daFile = open(path, 'r')
        headers = []
        self.das = []
        isFirst = True
        for l in daFile.readlines():
            da = [x for x in re.split('(?:\t+|\s{2,})', l) if len(x) > 0]
            # skip the empty line
            if len(da) < 3:
                continue

            if isFirst:
                headers = da
                isFirst = False
            else:
                daVal = {}
                # do converstions
                for i, h in enumerate(headers):
                    da[i] = da[i].strip()
                    if da[i] == 'yes':
                        da[i] = True
                    elif da[i] == 'no':
                        da[i] = False
                    else:
                        try:
                            da[i] = int(da[i])
                        except:
                            pass
                    daVal[h] = da[i]

                self.das.append(daVal)
        daFile.close()

    @staticmethod
    def instance():
        """Instance of donorAcceptor
        :return:
        """
        if DonAcceptor._instance is None:
            DonAcceptor._instance = DonAcceptor()
        return DonAcceptor._instance

    def __donors(self):
        return [a for a in self.das if a['Donor']]

    def __acceptors(self):
        return [a for a in self.das if a['Acceptor']]

    def showDonors(self):
        return [a for a in self.das if a['Donor']]

    @staticmethod
    def donors(atom, nextAtom):
        return DonAcceptor.instance().__probableDonors(atom, nextAtom)

    @staticmethod
    def acceptors(atom, nextAtom):
        return DonAcceptor.instance().__probableAcceptors(atom, nextAtom)

    def __probableDonors(self, atom, nextAtom):
        """ returns atom (and its orbital) and its hydrogen atom for donors"""
        for donor in self.__donors():
            if atom.symbol == donor['Atom'] and (donor['Residue'] == atom.residue or donor['Residue'] == 'Backbone'):
                hAtom = nextAtom(atom, 1)
                while hAtom.resId == atom.resId and hAtom.iCode == atom.iCode:
                    nAtom = nextAtom(hAtom, 1)
                    if nAtom == hAtom:
                        break
                    hAtom = nAtom
                    if hAtom.symbol == donor['Hydrogen']:
                        yield (atom, donor['orbital'], hAtom, donor)

    def __probableAcceptors(self, atom, nextAtom):
        for acceptor in self.__acceptors():
            if atom.symbol == acceptor['Atom'] and (
                            acceptor['Residue'] == atom.residue or acceptor['Residue'] == 'Backbone'):
                preAtom = nextAtom(atom, acceptor['preAtom(dist)'])
                yield (atom, acceptor['orbital'], preAtom, acceptor)
