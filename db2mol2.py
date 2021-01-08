#!/usr/bin/env python
import sys,os
import argparse

parser = argparse.ArgumentParser(description='Convert DB2 format compound library back to mol2.')
parser.add_argument('-i', '--input', type=str, required=True, help='Input db2 file. Could be compressed or uncompressed')
parser.add_argument('-o','--output', type=str, required=True, help='Output file')
parser.add_argument('-s','--outputsolv', type=str, help='Generate output.solv file', default='')

args = parser.parse_args()

infile = args.input
outfile = args.output
outputsolv = args.outputsolv

# Load all lines
if infile.endswith(".gz"):
    os.system("gunzip %s"%infile)
    infile = infile[:-3]

class Mblock:
    def __init__(self,lines):
        mlines = list()
        for line in lines:
            if line[0] == 'M':
                mlines.append(line)
        self.mlines = mlines

class Ablock:
    def __init__(self,lines):
        alines = list()
        for line in lines:
            if line[0] == 'A':
                alines.append(line)
        self.alines = alines 
        self.index = list()
        self.name = dict()
        self.atomtype = dict()
        self.vdwtype = dict()
        self.color = dict()
        self.properties = dict()
        self.charge = dict()
        for line in self.alines:
            items = line.split()
            index = int(items[1])
            self.index.append(index)
            self.name[index] = items[2]
            self.atomtype[index] = items[3]
            self.vdwtype[index] = int(items[4])
            self.color[index] = int(items[5])
            self.properties[index] = [ float(x) for x in items[6:] ]
            self.charge[index] = float(items[6])

class Bblock:
    def __init__(self,lines):
        blines = list()
        for line in lines:
            if line[0] == 'B':
                blines.append(line)
        self.blines = blines 

class Xblock:
    def __init__(self,lines):
        xlines = list()
        for line in lines:
            if line[0] == 'X':
                xlines.append(line)
        self.xlines = xlines
        # coordnum atomnum confnum x y z
        self.x = dict()
        for line in self.xlines:
            items = line.split()
            coordnum = int(items[1])
            atomnum = int(items[2])
            confnum = int(items[3])
            x = float(items[4])
            y = float(items[5])
            z = float(items[6])
            self.x[coordnum] = [atomnum,confnum,x,y,z]

class RBlock:
    def __init__(self,lines):
        rlines = list()
        for line in lines:
            if line[0] == 'R':
                rlines.append(line)
        self.rlines = rlines
        
class CBlock:
    def __init__(self,lines):
        clines = list()
        for line in lines:
            if line[0] == 'C':
                clines.append(line)
        self.clines = clines
        self.c = dict()
        for line in self.clines:
            items = line.split()
            # C confnum coordstart coordend
            confnum = int(items[1])
            coordstart = int(items[2])
            coordend = int(items[3])
            self.c[confnum] = [coordstart, coordend]

class SBlock:
    def __init__(self,lines):
        self.set = dict()
        self.slines = list()
        for line in lines:
            if line[0] == 'S':
                self.slines.append(line)

        for setblock in self.next_setblock():
            # print(">>>>>")
            # for line in setblock:
            #     print(line.strip())
            # print("<<<<<")
            # S setnum #lines #confs_total broken hydrogens omega_energy
            # S setnum linenum #confs confs [until full column]
            headline = setblock[0]
            otherlines = setblock[1:]
            items = headline.split()
            setnum = int(items[1])
            nlines = int(items[2])
            nconfs = int(items[3])
            iline = 0
            iconf = 0
            current_confs = list()
            for line in otherlines:
                iline += 1
                items = line.split()
                iconf += int(items[3])
                current_confs.extend( [ int(x) for x in items[4:] ] )
            assert iline == nlines
            assert iconf == nconfs == len(current_confs)
            self.set[setnum] = current_confs 
            # print(current_confs)

    def next_setblock(self):
        current_block = list()
        start = True
        linemax = 0
        linecount = 0
        for line in self.slines:
            if start:
                start = False
                current_block = list()
                current_block.append(line)
                items = line.split()
                linemax = int(items[2])
            else:
                linecount += 1
                if linecount != linemax:
                    current_block.append(line)
                else:
                    current_block.append(line)
                    start = True
                    linecount = 0
                    yield current_block

class DB2File:
    def __init__(self,infile):
        self.db2blocks = list()
        for lines in self.next_block(infile):
            # print(">>>>> Load one db2block")
            self.db2blocks.append(DB2Block(lines))

    def next_block(self,infile):
        block = list()
        for line in open(infile):
            if line[0] == 'E':
                yield block
                block = list()
            else:
                block.append(line)

    def convert_to_mol2(self,outfile):
        with open(outfile, 'w') as ofp:
            for db2block in self.db2blocks:
                db2block.write_mol2(ofp)

    def generate_outputsolv(self,outfile):
        outfilelines = list()
        headline = self._get_outputsolv_headline()
        atomline = self._get_outputsolv_atomline()
        with open(outfile, 'w') as ofp:
            ofp.write(headline)
            for tmp in atomline:
                ofp.write(tmp)

    def _get_outputsolv_headline(self):
        dbblock =self.db2blocks[0]
        m = dbblock.m
        atomnumber = int( m.mlines[0].split()[3] )
        secondline = m.mlines[1]
        items = secondline.split()
        charge = float(items[1])
        polsolv = float(items[2])
        apolsolv = float(items[3])
        totalsolv = float(items[4])
        surfa = float(items[5])
        mol2filename = 'output.mol2'
        output = f"{mol2filename} {atomnumber} {charge} {polsolv} {surfa} {apolsolv} {totalsolv}\n"
        return output

    def _get_outputsolv_atomline(self):
        dbblock =self.db2blocks[0]
        atomlines = dbblock.a.alines
        output = list()
        for line in atomlines:
            items = line.split()
            #A NUM NAME TYPEX DT CO +CHA.RGEX +POLAR.SOL +APOLA.SOL +TOTAL.SOL SURFA.REA
            charge = float(items[6])
            pol = float(items[7])
            apol = float(items[8])
            total = float(items[9])
            surfa = float(items[10])
            tmp = f"{charge} {pol} {total} {surfa} {apol}\n"
            output.append(tmp)
        return output

class DB2Block:
    def __init__(self,input_content,convert_back = True):
        if hasattr(input_content,'readlines'):
            self.lines = input_content.readlines()
        else:
            self.lines = input_content
        self.m = Mblock(self.lines)
        self.a = Ablock(self.lines)
        self.b = Bblock(self.lines)
        self.x = Xblock(self.lines)
        self.r = RBlock(self.lines)
        self.c = CBlock(self.lines)
        self.s = SBlock(self.lines)
        if convert_back:
            self.convert_mol2()

    def convert_mol2(self):
        self.mol2 = list()
        bondlines = self._get_bondlines()
        for setnum in self.s.set:
            thisset = self.s.set[setnum]
            atomlines = self._get_atomlines(thisset)
            headerlines = self._get_headerlines( annotations=None )
            mol2lines = self._get_mol2_lines(headerlines,atomlines,bondlines)
            self.mol2.append(mol2lines)

    def _get_atomlines(self, thisset):
        atomlines = list()
        allx = list()
        for confnum in thisset:
            conf = self.c.c[confnum]
            for xnum in range(conf[0],conf[1]+1):
                x = self.x.x[xnum]
                allx.append(x)
        allx = sorted(allx, key = lambda x:x[0])
        for x in allx:
            index = x[0]
            cx,cy,cz = x[2],x[3],x[4]
            name = self.a.name[index]
            atomtype = self.a.atomtype[index]
            resindex = 1 
            resname = "UNK"
            charge = self.a.charge[index]
            tmp = f"{index}\t{name}\t{cx}\t{cy}\t{cz}\t{atomtype}\t{resindex}\t{resname}\t{charge}\n"
            atomlines.append(tmp)
        return atomlines
        # [61, 69, 3.7622, -5.5219, -1.3007]

    def _get_bondlines(self):
        bondlines = list()
        for line in self.b.blines:
            bondlines.append(line[1:])
        return bondlines

    def _get_headerlines(self, annotations):
        headerlines = list()
        headerlines.append("5dby_DIF_1_ligand \n")
        headerlines.append("   29    30     1     0     0 \n")
        headerlines.append("SMALL \n")
        headerlines.append("NO_CHARGES \n")
        headerlines.append("\n")
        headerlines.append("\n")
        return headerlines

    def _get_mol2_lines(self,headerlines,atomlines,bondlines ):
        mol2lines = list()
        mol2lines.append("@<TRIPOS>MOLECULE\n")
        mol2lines.extend(headerlines)
        mol2lines.append("@<TRIPOS>ATOM\n")
        mol2lines.extend(atomlines)
        mol2lines.append("@<TRIPOS>BOND\n")
        mol2lines.extend(bondlines)
        return mol2lines

    def write_mol2(self,ofp):
        for mol2 in self.mol2:
            for line in mol2:
                ofp.write(line)

class Conformation:
    def __init__(self,lines):
        self.lines = lines
    def write_mol2(self,ofp):
        for line in self.lines:
            ofp.write_line(line)

a = DB2File(infile)
a.convert_to_mol2(outfile)
if outputsolv:
    a.generate_outputsolv(outputsolv)

