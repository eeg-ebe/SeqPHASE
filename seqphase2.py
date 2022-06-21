import sys

import math as python_lib_Math
import math as Math
import inspect as python_lib_Inspect
import sys as python_lib_Sys
import builtins as python_lib_Builtins
import functools as python_lib_Functools
import re as python_lib_Re
import traceback as python_lib_Traceback
from io import StringIO as python_lib_io_StringIO


class _hx_AnonObject:
    _hx_disable_getattr = False
    def __init__(self, fields):
        self.__dict__ = fields
    def __repr__(self):
        return repr(self.__dict__)
    def __contains__(self, item):
        return item in self.__dict__
    def __getitem__(self, item):
        return self.__dict__[item]
    def __getattr__(self, name):
        if (self._hx_disable_getattr):
            raise AttributeError('field does not exist')
        else:
            return None
    def _hx_hasattr(self,field):
        self._hx_disable_getattr = True
        try:
            getattr(self, field)
            self._hx_disable_getattr = False
            return True
        except AttributeError:
            self._hx_disable_getattr = False
            return False



class Enum:
    _hx_class_name = "Enum"
    __slots__ = ("tag", "index", "params")
    _hx_fields = ["tag", "index", "params"]
    _hx_methods = ["__str__"]

    def __init__(self,tag,index,params):
        self.tag = tag
        self.index = index
        self.params = params

    def __str__(self):
        if (self.params is None):
            return self.tag
        else:
            return self.tag + '(' + (', '.join(str(v) for v in self.params)) + ')'



class Class: pass


class EReg:
    _hx_class_name = "EReg"
    __slots__ = ("pattern", "matchObj", "_hx_global")
    _hx_fields = ["pattern", "matchObj", "global"]
    _hx_methods = ["split"]

    def __init__(self,r,opt):
        self.matchObj = None
        self._hx_global = False
        options = 0
        _g = 0
        _g1 = len(opt)
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            c = (-1 if ((i >= len(opt))) else ord(opt[i]))
            if (c == 109):
                options = (options | python_lib_Re.M)
            if (c == 105):
                options = (options | python_lib_Re.I)
            if (c == 115):
                options = (options | python_lib_Re.S)
            if (c == 117):
                options = (options | python_lib_Re.U)
            if (c == 103):
                self._hx_global = True
        self.pattern = python_lib_Re.compile(r,options)

    def split(self,s):
        if self._hx_global:
            ret = []
            lastEnd = 0
            x = python_HaxeIterator(python_lib_Re.finditer(self.pattern,s))
            while x.hasNext():
                x1 = x.next()
                x2 = HxString.substring(s,lastEnd,x1.start())
                ret.append(x2)
                lastEnd = x1.end()
            x = HxString.substr(s,lastEnd,None)
            ret.append(x)
            return ret
        else:
            self.matchObj = python_lib_Re.search(self.pattern,s)
            if (self.matchObj is None):
                return [s]
            else:
                return [HxString.substring(s,0,self.matchObj.start()), HxString.substr(s,self.matchObj.end(),None)]



class Individual:
    _hx_class_name = "Individual"
    __slots__ = ("name", "seq1", "seq2", "index", "p")
    _hx_fields = ["name", "seq1", "seq2", "index", "p"]
    _hx_methods = ["getName", "getSequence1", "getSequence2", "getIndex", "getProbability", "isHomozygous", "isHeterozygous", "getFasta"]

    def __init__(self,name,seq1,seq2,index,p):
        self.name = name
        self.seq1 = seq1
        self.seq2 = seq2
        self.index = index
        self.p = p

    def getName(self):
        return self.name

    def getSequence1(self):
        return self.seq1

    def getSequence2(self):
        return self.seq2

    def getIndex(self):
        return self.index

    def getProbability(self):
        return self.p

    def isHomozygous(self):
        return (self.seq1 == self.seq2)

    def isHeterozygous(self):
        return (self.seq1 != self.seq2)

    def getFasta(self,reduce,outFile):
        result = haxe_ds_List()
        if (self.p == 1.0):
            if (self.isHomozygous() and reduce):
                result.add((">" + HxOverrides.stringOrNull(self.name)))
                result.add(self.seq1)
            else:
                result.add(((">" + HxOverrides.stringOrNull(self.name)) + "a"))
                result.add(self.seq1)
                result.add(((">" + HxOverrides.stringOrNull(self.name)) + "b"))
                result.add(self.seq2)
        elif (self.isHomozygous() and reduce):
            result.add(((((((">" + HxOverrides.stringOrNull(self.name)) + "-") + Std.string(self.index)) + "(") + Std.string(self.p)) + ")"))
            result.add(self.seq1)
        else:
            result.add(((((((">" + HxOverrides.stringOrNull(self.name)) + "-") + Std.string(self.index)) + "a(") + Std.string(self.p)) + ")"))
            result.add(self.seq1)
            result.add(((((((">" + HxOverrides.stringOrNull(self.name)) + "-") + Std.string(self.index)) + "b(") + Std.string(self.p)) + ")"))
            result.add(self.seq2)
        return result.join("\n")



class Individuals:
    _hx_class_name = "Individuals"
    __slots__ = ("inds", "outFile")
    _hx_fields = ["inds", "outFile"]
    _hx_methods = ["addIndividual", "setOutFile", "getFasta"]

    def __init__(self):
        self.outFile = None
        self.inds = list()

    def addIndividual(self,ind):
        _this = self.inds
        _this.append(ind)

    def setOutFile(self,b):
        self.outFile = b

    def getFasta(self,reduce,sort):
        result = haxe_ds_List()
        def _hx_local_0(i1,i2):
            if (i1.getName() < i2.getName()):
                return -1
            if (i1.getName() > i2.getName()):
                return 1
            return (i1.getIndex() - i2.getIndex())
        self.inds.sort(key= python_lib_Functools.cmp_to_key(_hx_local_0))
        _g = 0
        _g1 = self.inds
        while (_g < len(_g1)):
            ind = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
            _g = (_g + 1)
            result.add(ind.getFasta(reduce,self.outFile))
        return result.join("\n")



class SeqPhase2:
    _hx_class_name = "SeqPhase2"
    __slots__ = ()
    _hx_statics = ["code", "supercode", "revSupercode", "removePossibleLineEndingR", "genConstString", "nSplit", "nrsToSeq", "parsePairFile", "stripBraces", "getSumCode", "processLines", "parseOutFile", "parse", "main"]

    @staticmethod
    def removePossibleLineEndingR(line):
        if (len(line) >= 1):
            if ((("" if ((0 >= len(line))) else line[0])) == "\r"):
                line = HxString.substr(line,1,None)
            index = (len(line) - 1)
            if ((("" if (((index < 0) or ((index >= len(line))))) else line[index])) == "\r"):
                line = HxString.substr(line,0,(len(line) - 1))
        return StringTools.trim(line)

    @staticmethod
    def genConstString(length):
        result = haxe_ds_List()
        _g = 0
        _g1 = length
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            result.add(".")
        return result.join("")

    @staticmethod
    def nSplit(s):
        result = haxe_ds_List()
        _g = 0
        _g1 = len(s)
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            char = ("" if (((i < 0) or ((i >= len(s))))) else s[i])
            if (char == " "):
                continue
            result.add(char)
        return result

    @staticmethod
    def nrsToSeq(nrs,const):
        nrs1 = SeqPhase2.nSplit(nrs)
        if ((const is None) or ((const == ""))):
            const = SeqPhase2.genConstString(nrs1.length)
        codeS = haxe_ds_List()
        _g_head = nrs1.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            part = val
            i = Std.parseInt(part)
            r = SeqPhase2.code.h.get(i,None)
            if (r is None):
                raise haxe_Exception.thrown((((("Error: The .out/.out_pairs file contains an allele not present in the input dataset (" + Std.string(i)) + ", ") + Std.string(nrs1)) + ")"))
            codeS.add(r)
        doesNotMatch = False
        result = haxe_ds_List()
        _g = 0
        _g1 = len(const)
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            _hx_chr = ("" if (((i < 0) or ((i >= len(const))))) else const[i])
            if (_hx_chr == "."):
                nextCodeS = codeS.pop()
                if (nextCodeS is None):
                    doesNotMatch = True
                    break
                result.add(nextCodeS)
            else:
                result.add(_hx_chr)
        if (doesNotMatch or ((codeS.length != 0))):
            raise haxe_Exception.thrown((((("Error: The data in the const and in the input file do not match; please check input data. (" + Std.string(doesNotMatch)) + ", ") + Std.string(codeS.length)) + ")"))
        return result.join("")

    @staticmethod
    def parsePairFile(fileContent,const = None):
        if (const is None):
            const = ""
        result = Individuals()
        result.setOutFile(False)
        lineNo = 0
        currentIndividualName = ""
        index = 0
        _g = 0
        _g1 = fileContent.split("\n")
        while (_g < len(_g1)):
            line = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
            _g = (_g + 1)
            line = SeqPhase2.removePossibleLineEndingR(line)
            lineNo = (lineNo + 1)
            if (line != ""):
                if (HxString.substr(line,0,5) == "IND: "):
                    currentIndividualName = HxString.substr(line,5,None)
                    index = 1
                else:
                    parts = line.split(" , ")
                    if (len(parts) != 3):
                        raise haxe_Exception.thrown((("Please check line " + Std.string(lineNo)) + " in the input file!"))
                    seq1 = SeqPhase2.nrsToSeq((parts[0] if 0 < len(parts) else None),const)
                    seq2 = SeqPhase2.nrsToSeq((parts[1] if 1 < len(parts) else None),const)
                    probability = Std.parseFloat((parts[2] if 2 < len(parts) else None))
                    ind = index
                    index = (index + 1)
                    ind1 = Individual(currentIndividualName,seq1,seq2,ind,probability)
                    result.addIndividual(ind1)
        return result

    @staticmethod
    def stripBraces(c):
        if ((("" if ((0 >= len(c))) else c[0])) == "("):
            c = HxString.substr(c,1,None)
        index = (len(c) - 1)
        if ((("" if (((index < 0) or ((index >= len(c))))) else c[index])) == ")"):
            c = HxString.substr(c,0,(len(c) - 1))
        return c

    @staticmethod
    def getSumCode(c1,c2):
        c1 = SeqPhase2.stripBraces(c1)
        c2 = SeqPhase2.stripBraces(c2)
        code1 = Std.parseInt(c1)
        code2 = Std.parseInt(c2)
        code1S = SeqPhase2.code.h.get(code1,None)
        code2S = SeqPhase2.code.h.get(code2,None)
        codeSum = (SeqPhase2.supercode.h.get(code1S,None) + SeqPhase2.supercode.h.get(code2S,None))
        return SeqPhase2.revSupercode.h.get(codeSum,None)

    @staticmethod
    def processLines(indName,seqLine1,seqLine2,const,lineNo):
        result = haxe_ds_List()
        r = EReg(" +","g")
        pLine1 = r.split(seqLine1)
        pLine2 = r.split(seqLine2)
        if (len(pLine1) != len(pLine2)):
            raise haxe_Exception.thrown((((("Error: Sequences have different lengths. Please check sequences of individual " + ("null" if indName is None else indName)) + " (line ") + Std.string(((lineNo - 2)))) + "ff.) in input file."))
        if ((const is None) or ((const == ""))):
            const = SeqPhase2.genConstString(len(pLine1))
        r1 = haxe_ds_List()
        r2 = haxe_ds_List()
        _g = 0
        _g1 = len(pLine1)
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            c1 = (pLine1[i] if i >= 0 and i < len(pLine1) else None)
            c2 = (pLine2[i] if i >= 0 and i < len(pLine2) else None)
            processed1 = False
            processed2 = False
            if ((("" if ((0 >= len(c1))) else c1[0])) == "["):
                r1.add("N")
                processed1 = True
            if ((("" if ((0 >= len(c2))) else c2[0])) == "["):
                r2.add("N")
                processed2 = True
            if ((("" if ((0 >= len(c1))) else c1[0])) == "("):
                toAdd = SeqPhase2.getSumCode(c1,c2)
                if ((toAdd is None) or ((toAdd == ""))):
                    raise haxe_Exception.thrown((((((((((("Unexpected input: Parsing error for " + ("null" if indName is None else indName)) + " ") + ("null" if toAdd is None else toAdd)) + " ") + Std.string(i)) + " ") + ("null" if c1 is None else c1)) + " ") + ("null" if c2 is None else c2)) + " (1)"))
                r1.add(toAdd)
                processed1 = True
            if ((("" if ((0 >= len(c2))) else c2[0])) == "("):
                toAdd1 = SeqPhase2.getSumCode(c2,c1)
                if ((toAdd1 is None) or ((toAdd1 == ""))):
                    raise haxe_Exception.thrown((((((((((("Unexpected input: Parsing error for " + ("null" if indName is None else indName)) + " ") + ("null" if toAdd1 is None else toAdd1)) + " ") + Std.string(i)) + " ") + ("null" if c1 is None else c1)) + " ") + ("null" if c2 is None else c2)) + " (2)"))
                r2.add(toAdd1)
                processed2 = True
            if (not processed1):
                _this = SeqPhase2.code
                key = Std.parseInt(c1)
                toAdd2 = _this.h.get(key,None)
                if ((toAdd2 is None) or ((toAdd2 == ""))):
                    raise haxe_Exception.thrown((((((((((("Unexpected input: Parsing error for " + ("null" if indName is None else indName)) + " ") + ("null" if toAdd2 is None else toAdd2)) + " ") + Std.string(i)) + " ") + ("null" if c1 is None else c1)) + " ") + ("null" if c2 is None else c2)) + " (3)"))
                r1.add(toAdd2)
            if (not processed2):
                _this1 = SeqPhase2.code
                key1 = Std.parseInt(c2)
                toAdd3 = _this1.h.get(key1,None)
                if ((toAdd3 is None) or ((toAdd3 == ""))):
                    raise haxe_Exception.thrown((((((((((("Unexpected input: Parsing error for " + ("null" if indName is None else indName)) + " ") + ("null" if toAdd3 is None else toAdd3)) + " ") + Std.string(i)) + " ") + ("null" if c1 is None else c1)) + " ") + ("null" if c2 is None else c2)) + " (4)"))
                r2.add(toAdd3)
        seq1 = haxe_ds_List()
        seq2 = haxe_ds_List()
        doesNotMatch = False
        _g = 0
        _g1 = len(const)
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            _hx_chr = ("" if (((i < 0) or ((i >= len(const))))) else const[i])
            if (_hx_chr == "."):
                r1S = r1.pop()
                r2S = r2.pop()
                if ((r1S is None) or ((r2S is None))):
                    doesNotMatch = True
                    break
                seq1.add(r1S)
                seq2.add(r2S)
            else:
                seq1.add(_hx_chr)
                seq2.add(_hx_chr)
        if ((doesNotMatch or ((r1.length != 0))) or ((r2.length != 0))):
            raise haxe_Exception.thrown((((((("Error: The data in the const and in the input file do not match; please check input data. (" + Std.string(doesNotMatch)) + ", ") + Std.string(r1.length)) + ", ") + Std.string(r2.length)) + ")"))
        ind = Individual(indName,seq1.join(""),seq2.join(""),-1,1.0)
        result.add(ind)
        return result

    @staticmethod
    def parseOutFile(fileContent,const = None):
        if (const is None):
            const = ""
        result = Individuals()
        result.setOutFile(True)
        lineNo = 0
        currentIndividualName = None
        currentSeq1 = None
        currentSeq2 = None
        index = 0
        seenBeginBestPairs1 = False
        _g = 0
        _g1 = fileContent.split("\n")
        while (_g < len(_g1)):
            line = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
            _g = (_g + 1)
            line = SeqPhase2.removePossibleLineEndingR(line)
            lineNo = (lineNo + 1)
            if ((line == "BEGIN BESTPAIRS1") or ((line == "BEGIN GENOTYPES"))):
                seenBeginBestPairs1 = True
            elif ((line == "END BESTPAIRS1") or ((line == "END GENOTYPES"))):
                break
            elif seenBeginBestPairs1:
                isThatAnEmptyLine = StringTools.trim(line)
                if (isThatAnEmptyLine == ""):
                    continue
                if (currentIndividualName is None):
                    currentIndividualName = HxOverrides.arrayGet(line.split(" "), 1)
                    if (currentIndividualName is None):
                        currentIndividualName = line
                elif (currentSeq1 is None):
                    currentSeq1 = line
                else:
                    currentSeq2 = line
                    _g_head = SeqPhase2.processLines(currentIndividualName,currentSeq1,currentSeq2,const,lineNo).h
                    while (_g_head is not None):
                        val = _g_head.item
                        _g_head = _g_head.next
                        ind = val
                        result.addIndividual(ind)
                    currentIndividualName = None
                    currentSeq1 = None
        return result

    @staticmethod
    def parse(fileContent,const = None):
        if (const is None):
            const = ""
        if (HxString.substr(fileContent,0,10) == "**********"):
            return SeqPhase2.parseOutFile(fileContent,const)
        elif (HxString.substr(fileContent,0,5) == "IND: "):
            return SeqPhase2.parsePairFile(fileContent,const)
        raise haxe_Exception.thrown("Input file not recognized (this is neither a .out nor a .out_pairs file generated by PHASE)")

    @staticmethod
    def main():
        Sys.stderr().writeString("SeqPHASE commmand-line version, Step 2: converting PHASE output files into FASTA alignments\n\n")
        Sys.stderr().writeString("Reference:\nFlot, J.-F. (2010) SeqPHASE: a web tool for interconverting PHASE input/output files and FASTA sequence alignments\nMolecular Ecology Ressources 10 (1): 162-166\n\n")
        Sys.stderr().writeString("Usage: perl seqphase2.pl -c <constant position file> -i <input=PHASE output file> -o <output=FASTA file name> [-r for reduced output] [-s for sorted output] \nThe constant position file (optional) was generated during Step 1 (without such file, the output FASTA will only contain the variable positions).\nInput=PHASE output file (compulsory) may be either .out or .out_pairs.\nOutput=FASTA file name (when not specified, default is 'phased.fasta'): if generated from an .out file, it will contain a list of phased haplotypes with 1-letter indetermination code letters (R, W, M, Y, S or K) at positions where phase certainty is inferior to a certain threshold (90% if PHASE default running options were used); if generated from an .out_pairs file, it will contain all possible haplotype pairs for each individual with their respective probabilities shown between parentheses.\nThe FASTA output may be reduced to show only posterior probabilities inferior to 1 and only one sequence per homozygote (-r switch) and/or be sorted alphabetically (-s switch).\n\n")
        myArgs = Sys.args()
        fileContent = None
        outfile = "phased.fasta"
        constFileContent = ""
        reduce = False
        sort = False
        i = 0
        while (i < len(myArgs)):
            current = (myArgs[i] if i >= 0 and i < len(myArgs) else None)
            if (current == "-r"):
                reduce = True
            elif (current == "-s"):
                sort = True
            elif (current == "-i"):
                i = (i + 1)
                myArgs1 = i
                fileContent = sys_io_File.getContent((myArgs[myArgs1] if myArgs1 >= 0 and myArgs1 < len(myArgs) else None))
            elif (current == "-o"):
                i = (i + 1)
                outfile1 = i
                outfile = (myArgs[outfile1] if outfile1 >= 0 and outfile1 < len(myArgs) else None)
            elif (current == "-c"):
                i = (i + 1)
                myArgs2 = i
                constFileContent = sys_io_File.getContent((myArgs[myArgs2] if myArgs2 >= 0 and myArgs2 < len(myArgs) else None))
            else:
                Sys.stderr().writeString((("Error: Unknown commandline option " + ("null" if current is None else current)) + "\n"))
                Sys.exit(1)
            i = (i + 1)
        if (fileContent is None):
            Sys.stderr().writeString("No input file specified.\n")
            Sys.exit(1)
        result = SeqPhase2.parse(fileContent,constFileContent)
        resultFileContent = result.getFasta(reduce,sort)
        fout = sys_io_File.write(outfile,False)
        fout.writeString(resultFileContent)
        fout.close()
        Sys.stderr().writeString("A FASTA alignments of phased haplotypes pairs")
        if sort:
            Sys.stderr().writeString(" (sorted alphabetically)")
        Sys.stderr().writeString(((" has been saved under " + ("null" if outfile is None else outfile)) + ".\n"))
        Sys.exit(0)


class Std:
    _hx_class_name = "Std"
    __slots__ = ()
    _hx_statics = ["isOfType", "string", "parseInt", "shortenPossibleNumber", "parseFloat"]

    @staticmethod
    def isOfType(v,t):
        if ((v is None) and ((t is None))):
            return False
        if (t is None):
            return False
        if ((type(t) == type) and (t == Dynamic)):
            return (v is not None)
        isBool = isinstance(v,bool)
        if (((type(t) == type) and (t == Bool)) and isBool):
            return True
        if ((((not isBool) and (not ((type(t) == type) and (t == Bool)))) and ((type(t) == type) and (t == Int))) and isinstance(v,int)):
            return True
        vIsFloat = isinstance(v,float)
        tmp = None
        tmp1 = None
        if (((not isBool) and vIsFloat) and ((type(t) == type) and (t == Int))):
            f = v
            tmp1 = (((f != Math.POSITIVE_INFINITY) and ((f != Math.NEGATIVE_INFINITY))) and (not python_lib_Math.isnan(f)))
        else:
            tmp1 = False
        if tmp1:
            tmp1 = None
            try:
                tmp1 = int(v)
            except BaseException as _g:
                None
                tmp1 = None
            tmp = (v == tmp1)
        else:
            tmp = False
        if ((tmp and ((v <= 2147483647))) and ((v >= -2147483648))):
            return True
        if (((not isBool) and ((type(t) == type) and (t == Float))) and isinstance(v,(float, int))):
            return True
        if ((type(t) == type) and (t == str)):
            return isinstance(v,str)
        isEnumType = ((type(t) == type) and (t == Enum))
        if ((isEnumType and python_lib_Inspect.isclass(v)) and hasattr(v,"_hx_constructs")):
            return True
        if isEnumType:
            return False
        isClassType = ((type(t) == type) and (t == Class))
        if ((((isClassType and (not isinstance(v,Enum))) and python_lib_Inspect.isclass(v)) and hasattr(v,"_hx_class_name")) and (not hasattr(v,"_hx_constructs"))):
            return True
        if isClassType:
            return False
        tmp = None
        try:
            tmp = isinstance(v,t)
        except BaseException as _g:
            None
            tmp = False
        if tmp:
            return True
        if python_lib_Inspect.isclass(t):
            cls = t
            loop = None
            def _hx_local_1(intf):
                f = (intf._hx_interfaces if (hasattr(intf,"_hx_interfaces")) else [])
                if (f is not None):
                    _g = 0
                    while (_g < len(f)):
                        i = (f[_g] if _g >= 0 and _g < len(f) else None)
                        _g = (_g + 1)
                        if (i == cls):
                            return True
                        else:
                            l = loop(i)
                            if l:
                                return True
                    return False
                else:
                    return False
            loop = _hx_local_1
            currentClass = v.__class__
            result = False
            while (currentClass is not None):
                if loop(currentClass):
                    result = True
                    break
                currentClass = python_Boot.getSuperClass(currentClass)
            return result
        else:
            return False

    @staticmethod
    def string(s):
        return python_Boot.toString1(s,"")

    @staticmethod
    def parseInt(x):
        if (x is None):
            return None
        try:
            return int(x)
        except BaseException as _g:
            None
            base = 10
            _hx_len = len(x)
            foundCount = 0
            sign = 0
            firstDigitIndex = 0
            lastDigitIndex = -1
            previous = 0
            _g = 0
            _g1 = _hx_len
            while (_g < _g1):
                i = _g
                _g = (_g + 1)
                c = (-1 if ((i >= len(x))) else ord(x[i]))
                if (((c > 8) and ((c < 14))) or ((c == 32))):
                    if (foundCount > 0):
                        return None
                    continue
                else:
                    c1 = c
                    if (c1 == 43):
                        if (foundCount == 0):
                            sign = 1
                        elif (not (((48 <= c) and ((c <= 57))))):
                            if (not (((base == 16) and ((((97 <= c) and ((c <= 122))) or (((65 <= c) and ((c <= 90))))))))):
                                break
                    elif (c1 == 45):
                        if (foundCount == 0):
                            sign = -1
                        elif (not (((48 <= c) and ((c <= 57))))):
                            if (not (((base == 16) and ((((97 <= c) and ((c <= 122))) or (((65 <= c) and ((c <= 90))))))))):
                                break
                    elif (c1 == 48):
                        if (not (((foundCount == 0) or (((foundCount == 1) and ((sign != 0))))))):
                            if (not (((48 <= c) and ((c <= 57))))):
                                if (not (((base == 16) and ((((97 <= c) and ((c <= 122))) or (((65 <= c) and ((c <= 90))))))))):
                                    break
                    elif ((c1 == 120) or ((c1 == 88))):
                        if ((previous == 48) and ((((foundCount == 1) and ((sign == 0))) or (((foundCount == 2) and ((sign != 0))))))):
                            base = 16
                        elif (not (((48 <= c) and ((c <= 57))))):
                            if (not (((base == 16) and ((((97 <= c) and ((c <= 122))) or (((65 <= c) and ((c <= 90))))))))):
                                break
                    elif (not (((48 <= c) and ((c <= 57))))):
                        if (not (((base == 16) and ((((97 <= c) and ((c <= 122))) or (((65 <= c) and ((c <= 90))))))))):
                            break
                if (((foundCount == 0) and ((sign == 0))) or (((foundCount == 1) and ((sign != 0))))):
                    firstDigitIndex = i
                foundCount = (foundCount + 1)
                lastDigitIndex = i
                previous = c
            if (firstDigitIndex <= lastDigitIndex):
                digits = HxString.substring(x,firstDigitIndex,(lastDigitIndex + 1))
                try:
                    return (((-1 if ((sign == -1)) else 1)) * int(digits,base))
                except BaseException as _g:
                    return None
            return None

    @staticmethod
    def shortenPossibleNumber(x):
        r = ""
        _g = 0
        _g1 = len(x)
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            c = ("" if (((i < 0) or ((i >= len(x))))) else x[i])
            _g2 = HxString.charCodeAt(c,0)
            if (_g2 is None):
                break
            else:
                _g3 = _g2
                if (((((((((((_g3 == 57) or ((_g3 == 56))) or ((_g3 == 55))) or ((_g3 == 54))) or ((_g3 == 53))) or ((_g3 == 52))) or ((_g3 == 51))) or ((_g3 == 50))) or ((_g3 == 49))) or ((_g3 == 48))) or ((_g3 == 46))):
                    r = (("null" if r is None else r) + ("null" if c is None else c))
                else:
                    break
        return r

    @staticmethod
    def parseFloat(x):
        try:
            return float(x)
        except BaseException as _g:
            None
            if (x is not None):
                r1 = Std.shortenPossibleNumber(x)
                if (r1 != x):
                    return Std.parseFloat(r1)
            return Math.NaN


class Float: pass


class Int: pass


class Bool: pass


class Dynamic: pass


class StringTools:
    _hx_class_name = "StringTools"
    __slots__ = ()
    _hx_statics = ["isSpace", "ltrim", "rtrim", "trim"]

    @staticmethod
    def isSpace(s,pos):
        if (((len(s) == 0) or ((pos < 0))) or ((pos >= len(s)))):
            return False
        c = HxString.charCodeAt(s,pos)
        if (not (((c > 8) and ((c < 14))))):
            return (c == 32)
        else:
            return True

    @staticmethod
    def ltrim(s):
        l = len(s)
        r = 0
        while ((r < l) and StringTools.isSpace(s,r)):
            r = (r + 1)
        if (r > 0):
            return HxString.substr(s,r,(l - r))
        else:
            return s

    @staticmethod
    def rtrim(s):
        l = len(s)
        r = 0
        while ((r < l) and StringTools.isSpace(s,((l - r) - 1))):
            r = (r + 1)
        if (r > 0):
            return HxString.substr(s,0,(l - r))
        else:
            return s

    @staticmethod
    def trim(s):
        return StringTools.ltrim(StringTools.rtrim(s))


class Sys:
    _hx_class_name = "Sys"
    __slots__ = ()
    _hx_statics = ["exit", "args", "stderr"]

    @staticmethod
    def exit(code):
        python_lib_Sys.exit(code)

    @staticmethod
    def args():
        argv = python_lib_Sys.argv
        return argv[1:None]

    @staticmethod
    def stderr():
        return python_io_IoTools.createFileOutputFromText(python_lib_Sys.stderr)


class haxe_IMap:
    _hx_class_name = "haxe.IMap"
    __slots__ = ()


class haxe_Exception(Exception):
    _hx_class_name = "haxe.Exception"
    __slots__ = ("_hx___nativeStack", "_hx___skipStack", "_hx___nativeException", "_hx___previousException")
    _hx_fields = ["__nativeStack", "__skipStack", "__nativeException", "__previousException"]
    _hx_methods = ["unwrap", "toString", "get_message", "get_native"]
    _hx_statics = ["caught", "thrown"]
    _hx_interfaces = []
    _hx_super = Exception


    def __init__(self,message,previous = None,native = None):
        self._hx___previousException = None
        self._hx___nativeException = None
        self._hx___nativeStack = None
        self._hx___skipStack = 0
        super().__init__(message)
        self._hx___previousException = previous
        if ((native is not None) and Std.isOfType(native,BaseException)):
            self._hx___nativeException = native
            self._hx___nativeStack = haxe_NativeStackTrace.exceptionStack()
        else:
            self._hx___nativeException = self
            infos = python_lib_Traceback.extract_stack()
            if (len(infos) != 0):
                infos.pop()
            infos.reverse()
            self._hx___nativeStack = infos

    def unwrap(self):
        return self._hx___nativeException

    def toString(self):
        return self.get_message()

    def get_message(self):
        return str(self)

    def get_native(self):
        return self._hx___nativeException

    @staticmethod
    def caught(value):
        if Std.isOfType(value,haxe_Exception):
            return value
        elif Std.isOfType(value,BaseException):
            return haxe_Exception(str(value),None,value)
        else:
            return haxe_ValueException(value,None,value)

    @staticmethod
    def thrown(value):
        if Std.isOfType(value,haxe_Exception):
            return value.get_native()
        elif Std.isOfType(value,BaseException):
            return value
        else:
            e = haxe_ValueException(value)
            e._hx___skipStack = (e._hx___skipStack + 1)
            return e



class haxe_NativeStackTrace:
    _hx_class_name = "haxe.NativeStackTrace"
    __slots__ = ()
    _hx_statics = ["saveStack", "exceptionStack"]

    @staticmethod
    def saveStack(exception):
        pass

    @staticmethod
    def exceptionStack():
        exc = python_lib_Sys.exc_info()
        if (exc[2] is not None):
            infos = python_lib_Traceback.extract_tb(exc[2])
            infos.reverse()
            return infos
        else:
            return []


class haxe_ValueException(haxe_Exception):
    _hx_class_name = "haxe.ValueException"
    __slots__ = ("value",)
    _hx_fields = ["value"]
    _hx_methods = ["unwrap"]
    _hx_statics = []
    _hx_interfaces = []
    _hx_super = haxe_Exception


    def __init__(self,value,previous = None,native = None):
        self.value = None
        super().__init__(Std.string(value),previous,native)
        self.value = value

    def unwrap(self):
        return self.value



class haxe_ds_IntMap:
    _hx_class_name = "haxe.ds.IntMap"
    __slots__ = ("h",)
    _hx_fields = ["h"]
    _hx_methods = ["set"]
    _hx_interfaces = [haxe_IMap]

    def __init__(self):
        self.h = dict()

    def set(self,key,value):
        self.h[key] = value



class haxe_ds_List:
    _hx_class_name = "haxe.ds.List"
    __slots__ = ("h", "q", "length")
    _hx_fields = ["h", "q", "length"]
    _hx_methods = ["add", "pop", "toString", "join"]

    def __init__(self):
        self.q = None
        self.h = None
        self.length = 0

    def add(self,item):
        x = haxe_ds__List_ListNode(item,None)
        if (self.h is None):
            self.h = x
        else:
            self.q.next = x
        self.q = x
        _hx_local_0 = self
        _hx_local_1 = _hx_local_0.length
        _hx_local_0.length = (_hx_local_1 + 1)
        _hx_local_1

    def pop(self):
        if (self.h is None):
            return None
        x = self.h.item
        self.h = self.h.next
        if (self.h is None):
            self.q = None
        _hx_local_0 = self
        _hx_local_1 = _hx_local_0.length
        _hx_local_0.length = (_hx_local_1 - 1)
        _hx_local_1
        return x

    def toString(self):
        s_b = python_lib_io_StringIO()
        first = True
        l = self.h
        s_b.write("{")
        while (l is not None):
            if first:
                first = False
            else:
                s_b.write(", ")
            s_b.write(Std.string(Std.string(l.item)))
            l = l.next
        s_b.write("}")
        return s_b.getvalue()

    def join(self,sep):
        s_b = python_lib_io_StringIO()
        first = True
        l = self.h
        while (l is not None):
            if first:
                first = False
            else:
                s_b.write(Std.string(sep))
            s_b.write(Std.string(l.item))
            l = l.next
        return s_b.getvalue()



class haxe_ds__List_ListNode:
    _hx_class_name = "haxe.ds._List.ListNode"
    __slots__ = ("item", "next")
    _hx_fields = ["item", "next"]

    def __init__(self,item,next):
        self.item = item
        self.next = next



class haxe_ds_StringMap:
    _hx_class_name = "haxe.ds.StringMap"
    __slots__ = ("h",)
    _hx_fields = ["h"]
    _hx_interfaces = [haxe_IMap]

    def __init__(self):
        self.h = dict()



class haxe_exceptions_PosException(haxe_Exception):
    _hx_class_name = "haxe.exceptions.PosException"
    __slots__ = ("posInfos",)
    _hx_fields = ["posInfos"]
    _hx_methods = ["toString"]
    _hx_statics = []
    _hx_interfaces = []
    _hx_super = haxe_Exception


    def __init__(self,message,previous = None,pos = None):
        self.posInfos = None
        super().__init__(message,previous)
        if (pos is None):
            self.posInfos = _hx_AnonObject({'fileName': "(unknown)", 'lineNumber': 0, 'className': "(unknown)", 'methodName': "(unknown)"})
        else:
            self.posInfos = pos

    def toString(self):
        return ((((((((("" + HxOverrides.stringOrNull(super().toString())) + " in ") + HxOverrides.stringOrNull(self.posInfos.className)) + ".") + HxOverrides.stringOrNull(self.posInfos.methodName)) + " at ") + HxOverrides.stringOrNull(self.posInfos.fileName)) + ":") + Std.string(self.posInfos.lineNumber))



class haxe_exceptions_NotImplementedException(haxe_exceptions_PosException):
    _hx_class_name = "haxe.exceptions.NotImplementedException"
    __slots__ = ()
    _hx_fields = []
    _hx_methods = []
    _hx_statics = []
    _hx_interfaces = []
    _hx_super = haxe_exceptions_PosException


    def __init__(self,message = None,previous = None,pos = None):
        if (message is None):
            message = "Not implemented"
        super().__init__(message,previous,pos)


class haxe_io_Bytes:
    _hx_class_name = "haxe.io.Bytes"
    __slots__ = ("length", "b")
    _hx_fields = ["length", "b"]
    _hx_statics = ["ofString"]

    def __init__(self,length,b):
        self.length = length
        self.b = b

    @staticmethod
    def ofString(s,encoding = None):
        b = bytearray(s,"UTF-8")
        return haxe_io_Bytes(len(b),b)


class haxe_io_Encoding(Enum):
    __slots__ = ()
    _hx_class_name = "haxe.io.Encoding"
    _hx_constructs = ["UTF8", "RawNative"]
haxe_io_Encoding.UTF8 = haxe_io_Encoding("UTF8", 0, ())
haxe_io_Encoding.RawNative = haxe_io_Encoding("RawNative", 1, ())

class haxe_io_Error(Enum):
    __slots__ = ()
    _hx_class_name = "haxe.io.Error"
    _hx_constructs = ["Blocked", "Overflow", "OutsideBounds", "Custom"]

    @staticmethod
    def Custom(e):
        return haxe_io_Error("Custom", 3, (e,))
haxe_io_Error.Blocked = haxe_io_Error("Blocked", 0, ())
haxe_io_Error.Overflow = haxe_io_Error("Overflow", 1, ())
haxe_io_Error.OutsideBounds = haxe_io_Error("OutsideBounds", 2, ())


class haxe_io_Output:
    _hx_class_name = "haxe.io.Output"
    __slots__ = ("bigEndian",)
    _hx_fields = ["bigEndian"]
    _hx_methods = ["writeByte", "writeBytes", "set_bigEndian", "writeFullBytes", "writeString"]

    def writeByte(self,c):
        raise haxe_exceptions_NotImplementedException(None,None,_hx_AnonObject({'fileName': "haxe/io/Output.hx", 'lineNumber': 47, 'className': "haxe.io.Output", 'methodName': "writeByte"}))

    def writeBytes(self,s,pos,_hx_len):
        if (((pos < 0) or ((_hx_len < 0))) or (((pos + _hx_len) > s.length))):
            raise haxe_Exception.thrown(haxe_io_Error.OutsideBounds)
        b = s.b
        k = _hx_len
        while (k > 0):
            self.writeByte(b[pos])
            pos = (pos + 1)
            k = (k - 1)
        return _hx_len

    def set_bigEndian(self,b):
        self.bigEndian = b
        return b

    def writeFullBytes(self,s,pos,_hx_len):
        while (_hx_len > 0):
            k = self.writeBytes(s,pos,_hx_len)
            pos = (pos + k)
            _hx_len = (_hx_len - k)

    def writeString(self,s,encoding = None):
        b = haxe_io_Bytes.ofString(s,encoding)
        self.writeFullBytes(b,0,b.length)



class haxe_iterators_ArrayIterator:
    _hx_class_name = "haxe.iterators.ArrayIterator"
    __slots__ = ("array", "current")
    _hx_fields = ["array", "current"]
    _hx_methods = ["hasNext", "next"]

    def __init__(self,array):
        self.current = 0
        self.array = array

    def hasNext(self):
        return (self.current < len(self.array))

    def next(self):
        def _hx_local_3():
            def _hx_local_2():
                _hx_local_0 = self
                _hx_local_1 = _hx_local_0.current
                _hx_local_0.current = (_hx_local_1 + 1)
                return _hx_local_1
            return python_internal_ArrayImpl._get(self.array, _hx_local_2())
        return _hx_local_3()



class python_Boot:
    _hx_class_name = "python.Boot"
    __slots__ = ()
    _hx_statics = ["keywords", "toString1", "fields", "simpleField", "getInstanceFields", "getSuperClass", "getClassFields", "prefixLength", "unhandleKeywords"]

    @staticmethod
    def toString1(o,s):
        if (o is None):
            return "null"
        if isinstance(o,str):
            return o
        if (s is None):
            s = ""
        if (len(s) >= 5):
            return "<...>"
        if isinstance(o,bool):
            if o:
                return "true"
            else:
                return "false"
        if (isinstance(o,int) and (not isinstance(o,bool))):
            return str(o)
        if isinstance(o,float):
            try:
                if (o == int(o)):
                    return str(Math.floor((o + 0.5)))
                else:
                    return str(o)
            except BaseException as _g:
                None
                return str(o)
        if isinstance(o,list):
            o1 = o
            l = len(o1)
            st = "["
            s = (("null" if s is None else s) + "\t")
            _g = 0
            _g1 = l
            while (_g < _g1):
                i = _g
                _g = (_g + 1)
                prefix = ""
                if (i > 0):
                    prefix = ","
                st = (("null" if st is None else st) + HxOverrides.stringOrNull(((("null" if prefix is None else prefix) + HxOverrides.stringOrNull(python_Boot.toString1((o1[i] if i >= 0 and i < len(o1) else None),s))))))
            st = (("null" if st is None else st) + "]")
            return st
        try:
            if hasattr(o,"toString"):
                return o.toString()
        except BaseException as _g:
            None
        if hasattr(o,"__class__"):
            if isinstance(o,_hx_AnonObject):
                toStr = None
                try:
                    fields = python_Boot.fields(o)
                    _g = []
                    _g1 = 0
                    while (_g1 < len(fields)):
                        f = (fields[_g1] if _g1 >= 0 and _g1 < len(fields) else None)
                        _g1 = (_g1 + 1)
                        x = ((("" + ("null" if f is None else f)) + " : ") + HxOverrides.stringOrNull(python_Boot.toString1(python_Boot.simpleField(o,f),(("null" if s is None else s) + "\t"))))
                        _g.append(x)
                    fieldsStr = _g
                    toStr = (("{ " + HxOverrides.stringOrNull(", ".join([x1 for x1 in fieldsStr]))) + " }")
                except BaseException as _g:
                    None
                    return "{ ... }"
                if (toStr is None):
                    return "{ ... }"
                else:
                    return toStr
            if isinstance(o,Enum):
                o1 = o
                l = len(o1.params)
                hasParams = (l > 0)
                if hasParams:
                    paramsStr = ""
                    _g = 0
                    _g1 = l
                    while (_g < _g1):
                        i = _g
                        _g = (_g + 1)
                        prefix = ""
                        if (i > 0):
                            prefix = ","
                        paramsStr = (("null" if paramsStr is None else paramsStr) + HxOverrides.stringOrNull(((("null" if prefix is None else prefix) + HxOverrides.stringOrNull(python_Boot.toString1(o1.params[i],s))))))
                    return (((HxOverrides.stringOrNull(o1.tag) + "(") + ("null" if paramsStr is None else paramsStr)) + ")")
                else:
                    return o1.tag
            if hasattr(o,"_hx_class_name"):
                if (o.__class__.__name__ != "type"):
                    fields = python_Boot.getInstanceFields(o)
                    _g = []
                    _g1 = 0
                    while (_g1 < len(fields)):
                        f = (fields[_g1] if _g1 >= 0 and _g1 < len(fields) else None)
                        _g1 = (_g1 + 1)
                        x = ((("" + ("null" if f is None else f)) + " : ") + HxOverrides.stringOrNull(python_Boot.toString1(python_Boot.simpleField(o,f),(("null" if s is None else s) + "\t"))))
                        _g.append(x)
                    fieldsStr = _g
                    toStr = (((HxOverrides.stringOrNull(o._hx_class_name) + "( ") + HxOverrides.stringOrNull(", ".join([x1 for x1 in fieldsStr]))) + " )")
                    return toStr
                else:
                    fields = python_Boot.getClassFields(o)
                    _g = []
                    _g1 = 0
                    while (_g1 < len(fields)):
                        f = (fields[_g1] if _g1 >= 0 and _g1 < len(fields) else None)
                        _g1 = (_g1 + 1)
                        x = ((("" + ("null" if f is None else f)) + " : ") + HxOverrides.stringOrNull(python_Boot.toString1(python_Boot.simpleField(o,f),(("null" if s is None else s) + "\t"))))
                        _g.append(x)
                    fieldsStr = _g
                    toStr = (((("#" + HxOverrides.stringOrNull(o._hx_class_name)) + "( ") + HxOverrides.stringOrNull(", ".join([x1 for x1 in fieldsStr]))) + " )")
                    return toStr
            if ((type(o) == type) and (o == str)):
                return "#String"
            if ((type(o) == type) and (o == list)):
                return "#Array"
            if callable(o):
                return "function"
            try:
                if hasattr(o,"__repr__"):
                    return o.__repr__()
            except BaseException as _g:
                None
            if hasattr(o,"__str__"):
                return o.__str__([])
            if hasattr(o,"__name__"):
                return o.__name__
            return "???"
        else:
            return str(o)

    @staticmethod
    def fields(o):
        a = []
        if (o is not None):
            if hasattr(o,"_hx_fields"):
                fields = o._hx_fields
                if (fields is not None):
                    return list(fields)
            if isinstance(o,_hx_AnonObject):
                d = o.__dict__
                keys = d.keys()
                handler = python_Boot.unhandleKeywords
                for k in keys:
                    if (k != '_hx_disable_getattr'):
                        a.append(handler(k))
            elif hasattr(o,"__dict__"):
                d = o.__dict__
                keys1 = d.keys()
                for k in keys1:
                    a.append(k)
        return a

    @staticmethod
    def simpleField(o,field):
        if (field is None):
            return None
        field1 = (("_hx_" + field) if ((field in python_Boot.keywords)) else (("_hx_" + field) if (((((len(field) > 2) and ((ord(field[0]) == 95))) and ((ord(field[1]) == 95))) and ((ord(field[(len(field) - 1)]) != 95)))) else field))
        if hasattr(o,field1):
            return getattr(o,field1)
        else:
            return None

    @staticmethod
    def getInstanceFields(c):
        f = (list(c._hx_fields) if (hasattr(c,"_hx_fields")) else [])
        if hasattr(c,"_hx_methods"):
            f = (f + c._hx_methods)
        sc = python_Boot.getSuperClass(c)
        if (sc is None):
            return f
        else:
            scArr = python_Boot.getInstanceFields(sc)
            scMap = set(scArr)
            _g = 0
            while (_g < len(f)):
                f1 = (f[_g] if _g >= 0 and _g < len(f) else None)
                _g = (_g + 1)
                if (not (f1 in scMap)):
                    scArr.append(f1)
            return scArr

    @staticmethod
    def getSuperClass(c):
        if (c is None):
            return None
        try:
            if hasattr(c,"_hx_super"):
                return c._hx_super
            return None
        except BaseException as _g:
            None
        return None

    @staticmethod
    def getClassFields(c):
        if hasattr(c,"_hx_statics"):
            x = c._hx_statics
            return list(x)
        else:
            return []

    @staticmethod
    def unhandleKeywords(name):
        if (HxString.substr(name,0,python_Boot.prefixLength) == "_hx_"):
            real = HxString.substr(name,python_Boot.prefixLength,None)
            if (real in python_Boot.keywords):
                return real
        return name


class python_HaxeIterator:
    _hx_class_name = "python.HaxeIterator"
    __slots__ = ("it", "x", "has", "checked")
    _hx_fields = ["it", "x", "has", "checked"]
    _hx_methods = ["next", "hasNext"]

    def __init__(self,it):
        self.checked = False
        self.has = False
        self.x = None
        self.it = it

    def next(self):
        if (not self.checked):
            self.hasNext()
        self.checked = False
        return self.x

    def hasNext(self):
        if (not self.checked):
            try:
                self.x = self.it.__next__()
                self.has = True
            except BaseException as _g:
                None
                if Std.isOfType(haxe_Exception.caught(_g).unwrap(),StopIteration):
                    self.has = False
                    self.x = None
                else:
                    raise _g
            self.checked = True
        return self.has



class python_internal_ArrayImpl:
    _hx_class_name = "python.internal.ArrayImpl"
    __slots__ = ()
    _hx_statics = ["_get"]

    @staticmethod
    def _get(x,idx):
        if ((idx > -1) and ((idx < len(x)))):
            return x[idx]
        else:
            return None


class HxOverrides:
    _hx_class_name = "HxOverrides"
    __slots__ = ()
    _hx_statics = ["eq", "stringOrNull", "arrayGet"]

    @staticmethod
    def eq(a,b):
        if (isinstance(a,list) or isinstance(b,list)):
            return a is b
        return (a == b)

    @staticmethod
    def stringOrNull(s):
        if (s is None):
            return "null"
        else:
            return s

    @staticmethod
    def arrayGet(a,i):
        if isinstance(a,list):
            x = a
            if ((i > -1) and ((i < len(x)))):
                return x[i]
            else:
                return None
        else:
            return a[i]


class python_internal_MethodClosure:
    _hx_class_name = "python.internal.MethodClosure"
    __slots__ = ("obj", "func")
    _hx_fields = ["obj", "func"]
    _hx_methods = ["__call__"]

    def __init__(self,obj,func):
        self.obj = obj
        self.func = func

    def __call__(self,*args):
        return self.func(self.obj,*args)



class HxString:
    _hx_class_name = "HxString"
    __slots__ = ()
    _hx_statics = ["charCodeAt", "substring", "substr"]

    @staticmethod
    def charCodeAt(s,index):
        if ((((s is None) or ((len(s) == 0))) or ((index < 0))) or ((index >= len(s)))):
            return None
        else:
            return ord(s[index])

    @staticmethod
    def substring(s,startIndex,endIndex = None):
        if (startIndex < 0):
            startIndex = 0
        if (endIndex is None):
            return s[startIndex:]
        else:
            if (endIndex < 0):
                endIndex = 0
            if (endIndex < startIndex):
                return s[endIndex:startIndex]
            else:
                return s[startIndex:endIndex]

    @staticmethod
    def substr(s,startIndex,_hx_len = None):
        if (_hx_len is None):
            return s[startIndex:]
        else:
            if (_hx_len == 0):
                return ""
            if (startIndex < 0):
                startIndex = (len(s) + startIndex)
                if (startIndex < 0):
                    startIndex = 0
            return s[startIndex:(startIndex + _hx_len)]


class python_io_NativeOutput(haxe_io_Output):
    _hx_class_name = "python.io.NativeOutput"
    __slots__ = ("stream",)
    _hx_fields = ["stream"]
    _hx_methods = ["close"]
    _hx_statics = []
    _hx_interfaces = []
    _hx_super = haxe_io_Output


    def __init__(self,stream):
        self.stream = None
        self.set_bigEndian(False)
        self.stream = stream
        if (not stream.writable()):
            raise haxe_Exception.thrown("Read only stream")

    def close(self):
        self.stream.close()



class python_io_NativeBytesOutput(python_io_NativeOutput):
    _hx_class_name = "python.io.NativeBytesOutput"
    __slots__ = ()
    _hx_fields = []
    _hx_methods = ["writeByte", "writeBytes"]
    _hx_statics = []
    _hx_interfaces = []
    _hx_super = python_io_NativeOutput


    def __init__(self,stream):
        super().__init__(stream)

    def writeByte(self,c):
        self.stream.write(bytearray([c]))

    def writeBytes(self,s,pos,_hx_len):
        return self.stream.write(s.b[pos:(pos + _hx_len)])



class python_io_IOutput:
    _hx_class_name = "python.io.IOutput"
    __slots__ = ("bigEndian",)
    _hx_fields = ["bigEndian"]
    _hx_methods = ["set_bigEndian", "writeByte", "writeBytes", "close", "writeFullBytes", "writeString"]


class python_io_IFileOutput:
    _hx_class_name = "python.io.IFileOutput"
    __slots__ = ()
    _hx_interfaces = [python_io_IOutput]


class python_io_FileBytesOutput(python_io_NativeBytesOutput):
    _hx_class_name = "python.io.FileBytesOutput"
    __slots__ = ()
    _hx_fields = []
    _hx_methods = []
    _hx_statics = []
    _hx_interfaces = [python_io_IFileOutput]
    _hx_super = python_io_NativeBytesOutput


    def __init__(self,stream):
        super().__init__(stream)


class python_io_NativeTextOutput(python_io_NativeOutput):
    _hx_class_name = "python.io.NativeTextOutput"
    __slots__ = ()
    _hx_fields = []
    _hx_methods = ["writeBytes", "writeByte"]
    _hx_statics = []
    _hx_interfaces = []
    _hx_super = python_io_NativeOutput


    def __init__(self,stream):
        super().__init__(stream)
        if (not stream.writable()):
            raise haxe_Exception.thrown("Read only stream")

    def writeBytes(self,s,pos,_hx_len):
        return self.stream.buffer.write(s.b[pos:(pos + _hx_len)])

    def writeByte(self,c):
        self.stream.write("".join(map(chr,[c])))



class python_io_FileTextOutput(python_io_NativeTextOutput):
    _hx_class_name = "python.io.FileTextOutput"
    __slots__ = ()
    _hx_fields = []
    _hx_methods = []
    _hx_statics = []
    _hx_interfaces = [python_io_IFileOutput]
    _hx_super = python_io_NativeTextOutput


    def __init__(self,stream):
        super().__init__(stream)


class python_io_IoTools:
    _hx_class_name = "python.io.IoTools"
    __slots__ = ()
    _hx_statics = ["createFileOutputFromText", "createFileOutputFromBytes"]

    @staticmethod
    def createFileOutputFromText(t):
        return sys_io_FileOutput(python_io_FileTextOutput(t))

    @staticmethod
    def createFileOutputFromBytes(t):
        return sys_io_FileOutput(python_io_FileBytesOutput(t))


class sys_io_File:
    _hx_class_name = "sys.io.File"
    __slots__ = ()
    _hx_statics = ["getContent", "write"]

    @staticmethod
    def getContent(path):
        f = python_lib_Builtins.open(path,"r",-1,"utf-8",None,"")
        content = f.read(-1)
        f.close()
        return content

    @staticmethod
    def write(path,binary = None):
        if (binary is None):
            binary = True
        mode = ("wb" if binary else "w")
        f = python_lib_Builtins.open(path,mode,-1,None,None,(None if binary else ""))
        if binary:
            return python_io_IoTools.createFileOutputFromBytes(f)
        else:
            return python_io_IoTools.createFileOutputFromText(f)


class sys_io_FileOutput(haxe_io_Output):
    _hx_class_name = "sys.io.FileOutput"
    __slots__ = ("impl",)
    _hx_fields = ["impl"]
    _hx_methods = ["set_bigEndian", "writeByte", "writeBytes", "close", "writeFullBytes", "writeString"]
    _hx_statics = []
    _hx_interfaces = []
    _hx_super = haxe_io_Output


    def __init__(self,impl):
        self.impl = impl

    def set_bigEndian(self,b):
        return self.impl.set_bigEndian(b)

    def writeByte(self,c):
        self.impl.writeByte(c)

    def writeBytes(self,s,pos,_hx_len):
        return self.impl.writeBytes(s,pos,_hx_len)

    def close(self):
        self.impl.close()

    def writeFullBytes(self,s,pos,_hx_len):
        self.impl.writeFullBytes(s,pos,_hx_len)

    def writeString(self,s,encoding = None):
        self.impl.writeString(s)


Math.NEGATIVE_INFINITY = float("-inf")
Math.POSITIVE_INFINITY = float("inf")
Math.NaN = float("nan")
Math.PI = python_lib_Math.pi

def _hx_init_SeqPhase2_code():
    def _hx_local_0():
        _g = haxe_ds_IntMap()
        _g.set(1,"A")
        _g.set(2,"C")
        _g.set(3,"G")
        _g.set(4,"T")
        _g.set(-1,"N")
        _g.set(0,"-")
        return _g
    return _hx_local_0()
SeqPhase2.code = _hx_init_SeqPhase2_code()
def _hx_init_SeqPhase2_supercode():
    def _hx_local_0():
        _g = haxe_ds_StringMap()
        _g.h["A"] = 1
        _g.h["T"] = 2
        _g.h["G"] = 4
        _g.h["C"] = 8
        _g.h["R"] = 5
        _g.h["Y"] = 10
        _g.h["M"] = 9
        _g.h["K"] = 6
        _g.h["W"] = 3
        _g.h["S"] = 12
        _g.h["V"] = 13
        _g.h["B"] = 14
        _g.h["H"] = 11
        _g.h["D"] = 7
        _g.h["N"] = 15
        return _g
    return _hx_local_0()
SeqPhase2.supercode = _hx_init_SeqPhase2_supercode()
def _hx_init_SeqPhase2_revSupercode():
    def _hx_local_0():
        _g = haxe_ds_IntMap()
        _g.set(1,"A")
        _g.set(2,"T")
        _g.set(4,"G")
        _g.set(8,"C")
        _g.set(5,"R")
        _g.set(10,"Y")
        _g.set(9,"M")
        _g.set(3,"W")
        _g.set(12,"S")
        _g.set(6,"K")
        _g.set(13,"V")
        _g.set(14,"B")
        _g.set(11,"H")
        _g.set(7,"D")
        _g.set(15,"N")
        return _g
    return _hx_local_0()
SeqPhase2.revSupercode = _hx_init_SeqPhase2_revSupercode()
python_Boot.keywords = set(["and", "del", "from", "not", "with", "as", "elif", "global", "or", "yield", "assert", "else", "if", "pass", "None", "break", "except", "import", "raise", "True", "class", "exec", "in", "return", "False", "continue", "finally", "is", "try", "def", "for", "lambda", "while"])
python_Boot.prefixLength = len("_hx_")

SeqPhase2.main()
