import sys

import math as python_lib_Math
import math as Math
import inspect as python_lib_Inspect
import sys as python_lib_Sys
import builtins as python_lib_Builtins
import functools as python_lib_Functools
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


class Entry:
    _hx_class_name = "Entry"
    __slots__ = ("name", "seq", "line")
    _hx_fields = ["name", "seq", "line"]
    _hx_methods = ["addToSeq", "getName", "getSeq", "getLineNo"]

    def __init__(self,line = None,name = None,seq = None):
        if (line is None):
            line = -1
        self.line = line
        self.name = name
        self.seq = seq

    def addToSeq(self,s):
        if ((self.seq is None) or ((self.seq == ""))):
            self.seq = s
        else:
            _hx_local_0 = self
            _hx_local_1 = _hx_local_0.seq
            _hx_local_0.seq = (("null" if _hx_local_1 is None else _hx_local_1) + ("null" if s is None else s))
            _hx_local_0.seq

    def getName(self):
        return self.name

    def getSeq(self):
        if (self.seq is None):
            return ""
        return self.seq

    def getLineNo(self):
        return self.line



class FastaAlignmentParser:
    _hx_class_name = "FastaAlignmentParser"
    __slots__ = ("fastaContent", "seqLength")
    _hx_fields = ["fastaContent", "seqLength"]
    _hx_methods = ["getSeqLength", "getSequences"]
    _hx_statics = ["authorizedCharacters", "startsWith", "isWhitespace", "stripStringBegin", "stripStringEnd", "stripString"]

    def __init__(self,fileContent,allChecks,allSort,fileNr):
        self.fastaContent = list()
        self.seqLength = -1
        if (fileContent is None):
            return
        end = len(fileContent)
        while True:
            tmp = None
            if (end > 0):
                cCode = HxString.charCodeAt(fileContent,(end - 1))
                result = False
                _g = 0
                _g1 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                while (_g < len(_g1)):
                    ele = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                    _g = (_g + 1)
                    if (ele == cCode):
                        result = True
                        break
                tmp = result
            else:
                tmp = False
            if (not tmp):
                break
            end = (end - 1)
        s = HxString.substring(fileContent,0,end)
        begin = 0
        sLen = len(s)
        while True:
            tmp = None
            if (begin < sLen):
                cCode = HxString.charCodeAt(s,begin)
                result = False
                _g = 0
                _g1 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                while (_g < len(_g1)):
                    ele = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                    _g = (_g + 1)
                    if (ele == cCode):
                        result = True
                        break
                tmp = result
            else:
                tmp = False
            if (not tmp):
                break
            begin = (begin + 1)
        if (HxString.substr(s,begin,None) == ""):
            SeqPhase1Result.instance().addErr("Empty file!",fileNr)
            return
        if ((HxString.substr(fileContent,0,len(">")) != ">") and ((HxString.substr(fileContent,0,len(";")) != ";"))):
            SeqPhase1Result.instance().addErr("File does not seem to be a fasta file!",fileNr)
            return
        lines = fileContent.split("\n")
        if (len(lines) == 0):
            SeqPhase1Result.instance().addErr("Not a fasta but an empty file!",fileNr)
            return
        elif (len(lines) == 1):
            SeqPhase1Result.instance().addErr("Only 1 line detected! Please check data format (opening the alignment in MEGA (http://www.megasoftware.net/) and exporting it as FASTA again may solve the problem; alternatively, there may be a problem with end-of-line characters - see http://en.wikipedia.org/wiki/Newline for details!",fileNr)
            return
        entryMap = haxe_ds_StringMap()
        current = None
        lineNo = 0
        underscoreWarningOutputted = False
        _g = 0
        while (_g < len(lines)):
            line = (lines[_g] if _g >= 0 and _g < len(lines) else None)
            _g = (_g + 1)
            lineNo = (lineNo + 1)
            end = len(line)
            while True:
                line1 = None
                if (end > 0):
                    cCode = HxString.charCodeAt(line,(end - 1))
                    result = False
                    _g1 = 0
                    _g2 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                    while (_g1 < len(_g2)):
                        ele = (_g2[_g1] if _g1 >= 0 and _g1 < len(_g2) else None)
                        _g1 = (_g1 + 1)
                        if (ele == cCode):
                            result = True
                            break
                    line1 = result
                else:
                    line1 = False
                if (not line1):
                    break
                end = (end - 1)
            s = HxString.substring(line,0,end)
            begin = 0
            sLen = len(s)
            while True:
                line2 = None
                if (begin < sLen):
                    cCode1 = HxString.charCodeAt(s,begin)
                    result1 = False
                    _g3 = 0
                    _g4 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                    while (_g3 < len(_g4)):
                        ele1 = (_g4[_g3] if _g3 >= 0 and _g3 < len(_g4) else None)
                        _g3 = (_g3 + 1)
                        if (ele1 == cCode1):
                            result1 = True
                            break
                    line2 = result1
                else:
                    line2 = False
                if (not line2):
                    break
                begin = (begin + 1)
            line = HxString.substr(s,begin,None)
            if (HxString.substr(line,0,len(";")) == ";"):
                continue
            elif (HxString.substr(line,0,len(">")) == ">"):
                s1 = HxString.substr(line,1,None)
                end1 = len(s1)
                while True:
                    tmp = None
                    if (end1 > 0):
                        cCode2 = HxString.charCodeAt(s1,(end1 - 1))
                        result2 = False
                        _g5 = 0
                        _g6 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                        while (_g5 < len(_g6)):
                            ele2 = (_g6[_g5] if _g5 >= 0 and _g5 < len(_g6) else None)
                            _g5 = (_g5 + 1)
                            if (ele2 == cCode2):
                                result2 = True
                                break
                        tmp = result2
                    else:
                        tmp = False
                    if (not tmp):
                        break
                    end1 = (end1 - 1)
                s2 = HxString.substring(s1,0,end1)
                begin1 = 0
                sLen1 = len(s2)
                while True:
                    tmp1 = None
                    if (begin1 < sLen1):
                        cCode3 = HxString.charCodeAt(s2,begin1)
                        result3 = False
                        _g7 = 0
                        _g8 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                        while (_g7 < len(_g8)):
                            ele3 = (_g8[_g7] if _g7 >= 0 and _g7 < len(_g8) else None)
                            _g7 = (_g7 + 1)
                            if (ele3 == cCode3):
                                result3 = True
                                break
                        tmp1 = result3
                    else:
                        tmp1 = False
                    if (not tmp1):
                        break
                    begin1 = (begin1 + 1)
                indName = HxString.substr(s2,begin1,None)
                indNameCor = StringTools.replace(indName," ","_")
                if (indName != indNameCor):
                    if (not underscoreWarningOutputted):
                        SeqPhase1Result.instance().addWrn("Warning: PHASE does not accept spaces in individual names. These spaces got replaced by underscore characters.",fileNr)
                        underscoreWarningOutputted = True
                    indName = indNameCor
                if (indName in entryMap.h):
                    SeqPhase1Result.instance().addErr((((((("Repeat of name " + ("null" if indName is None else indName)) + " encountered in alignment (line ") + Std.string(lineNo)) + ", line ") + Std.string(entryMap.h.get(indName,None).getLineNo())) + ")"),fileNr)
                current = Entry(lineNo,indName)
                entryMap.h[indName] = current
                if (len(indName) == 0):
                    SeqPhase1Result.instance().addErr(("Missing sequence name, line " + Std.string(lineNo)),fileNr)
            else:
                _g9 = 0
                _g10 = len(line)
                while (_g9 < _g10):
                    i = _g9
                    _g9 = (_g9 + 1)
                    char = ("" if (((i < 0) or ((i >= len(line))))) else line[i]).upper()
                    if (not (char in FastaAlignmentParser.authorizedCharacters.h)):
                        SeqPhase1Result.instance().addErr(((((("Unknown character state " + ("null" if char is None else char)) + " in ") + HxOverrides.stringOrNull(current.getName())) + ", line ") + Std.string(lineNo)),fileNr)
                    elif (allChecks and ((FastaAlignmentParser.authorizedCharacters.h.get(char,None) == False))):
                        SeqPhase1Result.instance().addErr(((((("Unallowed state " + ("null" if char is None else char)) + " in ") + HxOverrides.stringOrNull(current.getName())) + ", line ") + Std.string(lineNo)),fileNr)
                _this = line.split("?")
                line = "N".join([python_Boot.toString1(x1,'') for x1 in _this])
                current.addToSeq(line.upper())
        if (current is None):
            SeqPhase1Result.instance().addErr("Corrupted Fasta File",fileNr)
            return
        if (current.getSeq() == ""):
            SeqPhase1Result.instance().addErr(((("Empty sequence " + HxOverrides.stringOrNull(current.getName())) + ", line ") + Std.string(current.getLineNo())),fileNr)
            return
        self.seqLength = len(current.getSeq())
        key = entryMap.keys()
        while key.hasNext():
            key1 = key.next()
            val = entryMap.h.get(key1,None)
            _this = self.fastaContent
            _this.append(val)
            if (len(val.getSeq()) != len(current.getSeq())):
                SeqPhase1Result.instance().addErr(((((((((((("Not all sequences in this file have equal lengths. E.g. sequence " + HxOverrides.stringOrNull(val.getName())) + " (line ") + Std.string(val.getLineNo())) + ") is of length ") + Std.string(len(val.getSeq()))) + " while sequence ") + HxOverrides.stringOrNull(current.getName())) + " (line ") + Std.string(current.getLineNo())) + ") is of length ") + Std.string(len(current.getSeq()))),fileNr)
                return
            if (HxString.substr(val.getSeq(),0,len("-")) == "-"):
                SeqPhase1Result.instance().addWrn((((("Sequence " + HxOverrides.stringOrNull(val.getName())) + " (line ") + Std.string(val.getLineNo())) + ") starts with '-'. Is it a real indel or did you mean 'N' or '?' (missing data)?"),fileNr)
            _this1 = val.getSeq()
            index = (len(val.getSeq()) - 1)
            if ((("" if (((index < 0) or ((index >= len(_this1))))) else _this1[index])) == "-"):
                SeqPhase1Result.instance().addWrn((((("Sequence " + HxOverrides.stringOrNull(val.getName())) + " (line ") + Std.string(val.getLineNo())) + ") ends with '-'. Is it a real indel or did you mean 'N' or '?' (missing data)?"),fileNr)
        def _hx_local_14(e1,e2):
            if allSort:
                nameE1 = e1.getName()
                nameE2 = e1.getName()
                aNameE1 = HxString.substring(nameE1,0,(len(nameE1) - 1))
                aNameE2 = HxString.substring(nameE2,0,(len(nameE2) - 1))
                if (aNameE1 < aNameE2):
                    return -1
                if (aNameE1 > aNameE2):
                    return 1
                index = (len(nameE1) - 1)
                lNameE1 = ("" if (((index < 0) or ((index >= len(nameE1))))) else nameE1[index])
                index = (len(nameE2) - 1)
                lNameE2 = ("" if (((index < 0) or ((index >= len(nameE2))))) else nameE2[index])
                if (lNameE1 < lNameE2):
                    return -1
                if (lNameE1 > lNameE2):
                    return 1
                return 0
            else:
                if (e1.getName() < e2.getName()):
                    return -1
                if (e1.getName() > e2.getName()):
                    return 1
                return 0
        self.fastaContent.sort(key= python_lib_Functools.cmp_to_key(_hx_local_14))
        if allChecks:
            if (HxOverrides.mod(len(self.fastaContent), 2) != 0):
                SeqPhase1Result.instance().addErr("Uneven number of sequences in alignment: please check data.",fileNr)
            lastName = None
            _g = 0
            _g1 = self.fastaContent
            while (_g < len(_g1)):
                entry = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                _g = (_g + 1)
                if (lastName is None):
                    lastName = entry.getName()
                    lastName = HxString.substring(lastName,0,(len(lastName) - 1))
                else:
                    curName = entry.getName()
                    curName = HxString.substring(curName,0,(len(curName) - 1))
                    if (lastName == curName):
                        lastName = None
            _hx_map = haxe_ds_StringMap()
            _g = 0
            _g1 = self.fastaContent
            while (_g < len(_g1)):
                entry = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                _g = (_g + 1)
                name = entry.getName()
                name = HxString.substring(name,0,(len(name) - 1))
                if (name in _hx_map.h):
                    i = _hx_map.h.get(name,None)
                    _hx_map.h[name] = (i + 1)
                else:
                    _hx_map.h[name] = 1
            name = _hx_map.keys()
            while name.hasNext():
                name1 = name.next()
                i = _hx_map.h.get(name1,None)
                if (i != 2):
                    lst = haxe_ds_List()
                    _g = 0
                    _g1 = self.fastaContent
                    while (_g < len(_g1)):
                        entry = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                        _g = (_g + 1)
                        nameEntry = entry.getName()
                        subnameEntry = HxString.substring(nameEntry,0,(len(nameEntry) - 1))
                        if (subnameEntry == name1):
                            lst.add(nameEntry)
                    if (i == 1):
                        SeqPhase1Result.instance().addErr((((("Found only one haplotype sequence (" + HxOverrides.stringOrNull(lst.join(","))) + ") for individual '") + ("null" if name1 is None else name1)) + "'!"),fileNr)
                    else:
                        SeqPhase1Result.instance().addErr((((((("Found " + Std.string(i)) + " haplotype sequences (") + HxOverrides.stringOrNull(lst.join(","))) + ") for individual '") + ("null" if name1 is None else name1)) + "'!"),fileNr)

    def getSeqLength(self):
        return self.seqLength

    def getSequences(self):
        return self.fastaContent

    @staticmethod
    def startsWith(t,s):
        return (HxString.substr(t,0,len(s)) == s)

    @staticmethod
    def isWhitespace(s,pos):
        cCode = HxString.charCodeAt(s,pos)
        result = False
        _g = 0
        _g1 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
        while (_g < len(_g1)):
            ele = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
            _g = (_g + 1)
            if (ele == cCode):
                result = True
                break
        return result

    @staticmethod
    def stripStringBegin(s):
        begin = 0
        sLen = len(s)
        while True:
            tmp = None
            if (begin < sLen):
                cCode = HxString.charCodeAt(s,begin)
                result = False
                _g = 0
                _g1 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                while (_g < len(_g1)):
                    ele = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                    _g = (_g + 1)
                    if (ele == cCode):
                        result = True
                        break
                tmp = result
            else:
                tmp = False
            if (not tmp):
                break
            begin = (begin + 1)
        return HxString.substr(s,begin,None)

    @staticmethod
    def stripStringEnd(s):
        end = len(s)
        while True:
            tmp = None
            if (end > 0):
                cCode = HxString.charCodeAt(s,(end - 1))
                result = False
                _g = 0
                _g1 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                while (_g < len(_g1)):
                    ele = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                    _g = (_g + 1)
                    if (ele == cCode):
                        result = True
                        break
                tmp = result
            else:
                tmp = False
            if (not tmp):
                break
            end = (end - 1)
        return HxString.substring(s,0,end)

    @staticmethod
    def stripString(s):
        end = len(s)
        while True:
            tmp = None
            if (end > 0):
                cCode = HxString.charCodeAt(s,(end - 1))
                result = False
                _g = 0
                _g1 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                while (_g < len(_g1)):
                    ele = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                    _g = (_g + 1)
                    if (ele == cCode):
                        result = True
                        break
                tmp = result
            else:
                tmp = False
            if (not tmp):
                break
            end = (end - 1)
        s1 = HxString.substring(s,0,end)
        begin = 0
        sLen = len(s1)
        while True:
            tmp = None
            if (begin < sLen):
                cCode = HxString.charCodeAt(s1,begin)
                result = False
                _g = 0
                _g1 = [9, 10, 11, 12, 13, 32, 133, 160, 5760, 8192, 8192, 8193, 8194, 8195, 8196, 8197, 8198, 8199, 8200, 8201, 8202, 8232, 8233, 8239, 8287, 12288, 6158, 8203, 8204, 8205, 8288, 65279]
                while (_g < len(_g1)):
                    ele = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
                    _g = (_g + 1)
                    if (ele == cCode):
                        result = True
                        break
                tmp = result
            else:
                tmp = False
            if (not tmp):
                break
            begin = (begin + 1)
        return HxString.substr(s1,begin,None)



class SeqPhase1:
    _hx_class_name = "SeqPhase1"
    __slots__ = ()
    _hx_statics = ["map1", "map2", "code", "doIt", "makeStr", "main"]

    @staticmethod
    def doIt(align1,align2,align3):
        SeqPhase1Result.instance().clear()
        al1 = FastaAlignmentParser(align1,False,False,1)
        al2 = FastaAlignmentParser(align2,True,True,2)
        al3 = FastaAlignmentParser(align3,True,False,3)
        expectedLength = (al1.getSeqLength() if ((al1.getSeqLength() > al2.getSeqLength())) else al2.getSeqLength())
        if (expectedLength <= al3.getSeqLength()):
            expectedLength = al3.getSeqLength()
        diffLength = False
        if ((al1.getSeqLength() != -1) and ((al1.getSeqLength() != expectedLength))):
            diffLength = True
        if ((al2.getSeqLength() != -1) and ((al2.getSeqLength() != expectedLength))):
            diffLength = True
        if ((al3.getSeqLength() != -1) and ((al3.getSeqLength() != expectedLength))):
            diffLength = True
        if diffLength:
            SeqPhase1Result.instance().addGeneralError("Not all input sequences have equal lengths, please check whether this is expected.")
        elif ((expectedLength == -1) or ((expectedLength == 0))):
            SeqPhase1Result.instance().addGeneralError("It seems that all given sequences are empty ...")
        if SeqPhase1Result.instance().hasErrors():
            return SeqPhase1Result.instance()
        al1a = haxe_ds_List()
        al1b = haxe_ds_List()
        _g = 0
        _g1 = al1.getSequences()
        while (_g < len(_g1)):
            entry = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
            _g = (_g + 1)
            seq1a = haxe_ds_List()
            seq1b = haxe_ds_List()
            _g2 = 0
            _g3 = len(entry.getSeq())
            while (_g2 < _g3):
                i = _g2
                _g2 = (_g2 + 1)
                _this = entry.getSeq()
                c = ("" if (((i < 0) or ((i >= len(_this))))) else _this[i])
                if (c in SeqPhase1.map1.h):
                    seq1a.add(SeqPhase1.map1.h.get(c,None))
                    seq1b.add(SeqPhase1.map2.h.get(c,None))
                else:
                    seq1a.add(c)
                    seq1b.add(c)
            al1a.add(Entry(entry.getLineNo(),entry.getName(),seq1a.join("")))
            al1b.add(Entry(entry.getLineNo(),entry.getName(),seq1b.join("")))
        varpos = haxe_ds_List()
        multipos = haxe_ds_List()
        multiposMap = haxe_ds_IntMap()
        constFileContent = haxe_ds_List()
        _g = 0
        _g1 = expectedLength
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            m = haxe_ds_StringMap()
            _g2_head = al1a.h
            while (_g2_head is not None):
                val = _g2_head.item
                _g2_head = _g2_head.next
                entry = val
                _this = entry.getSeq()
                key = ("" if (((i < 0) or ((i >= len(_this))))) else _this[i])
                m.h[key] = False
            _g3_head = al1b.h
            while (_g3_head is not None):
                val1 = _g3_head.item
                _g3_head = _g3_head.next
                entry1 = val1
                _this1 = entry1.getSeq()
                key1 = ("" if (((i < 0) or ((i >= len(_this1))))) else _this1[i])
                m.h[key1] = False
            _g2 = 0
            _g3 = al2.getSequences()
            while (_g2 < len(_g3)):
                entry2 = (_g3[_g2] if _g2 >= 0 and _g2 < len(_g3) else None)
                _g2 = (_g2 + 1)
                _this2 = entry2.getSeq()
                key2 = ("" if (((i < 0) or ((i >= len(_this2))))) else _this2[i])
                m.h[key2] = False
            _g4 = 0
            _g5 = al3.getSequences()
            while (_g4 < len(_g5)):
                entry3 = (_g5[_g4] if _g4 >= 0 and _g4 < len(_g5) else None)
                _g4 = (_g4 + 1)
                _this3 = entry3.getSeq()
                key3 = ("" if (((i < 0) or ((i >= len(_this3))))) else _this3[i])
                m.h[key3] = False
            mapLen = 0
            mapLenWithoutNs = 0
            lastKey = None
            key4 = m.keys()
            while key4.hasNext():
                key5 = key4.next()
                lastKey = key5
                mapLen = (mapLen + 1)
                if (key5 != "N"):
                    mapLenWithoutNs = (mapLenWithoutNs + 1)
            if (mapLen == 0):
                SeqPhase1Result.instance().addGeneralError("Bug detected: There seems to be a bug in this implementation of SeqPHASE. Please contact the author so that they can fix this bug.")
                constFileContent.add("X")
            elif (mapLen == 1):
                if (lastKey == "N"):
                    SeqPhase1Result.instance().addGeneralWarn((("Found only N/? s at position " + Std.string(((i + 1)))) + "."))
                elif (lastKey == "-"):
                    SeqPhase1Result.instance().addGeneralWarn((("Found only -'s at position " + Std.string(((i + 1)))) + ". This may indicate an alignment problem. Consider to recheck your input data!"))
                constFileContent.add(lastKey)
            else:
                constFileContent.add(".")
                varpos.add(i)
                if (mapLenWithoutNs > 2):
                    multipos.add(i)
                    multiposMap.set(i,True)
                else:
                    multiposMap.set(i,False)
        if (varpos.length == 0):
            SeqPhase1Result.instance().addGeneralError("Not a single variable position detected in dataset! Please check data.")
        else:
            SeqPhase1Result.instance().addNote((((("There are " + Std.string(varpos.length)) + " variable positions in your dataset, including ") + Std.string(multipos.length)) + " position(s) with more than two different states."))
        SeqPhase1Result.instance().setConstFile(constFileContent.join(""))
        lines = haxe_ds_List()
        x = ((al1a.length + ((len(al2.getSequences()) / 2))) + ((len(al3.getSequences()) / 2)))
        tmp = None
        try:
            tmp = int(x)
        except BaseException as _g:
            None
            tmp = None
        lines.add(("" + Std.string(tmp)))
        lines.add(("" + Std.string(varpos.length)))
        l1 = haxe_ds_List()
        l2 = haxe_ds_List()
        _g4_head = varpos.h
        while (_g4_head is not None):
            val = _g4_head.item
            _g4_head = _g4_head.next
            pos = val
            l1.add(("" + Std.string(((pos + 1)))))
            if multiposMap.h.get(pos,None):
                l2.add("M")
            else:
                l2.add("S")
        lines.add(("P " + HxOverrides.stringOrNull(l1.join(" "))))
        lines.add((HxOverrides.stringOrNull(l2.join(" ")) + " "))
        it1_head = al1a.h
        it2_head = al1b.h
        while (it1_head is not None):
            val = it1_head.item
            it1_head = it1_head.next
            e1 = val
            val1 = it2_head.item
            it2_head = it2_head.next
            e2 = val1
            lines.add(e1.getName())
            line1 = haxe_ds_List()
            line2 = haxe_ds_List()
            _g5_head = varpos.h
            while (_g5_head is not None):
                val2 = _g5_head.item
                _g5_head = _g5_head.next
                i = val2
                _this = e1.getSeq()
                char = ("" if (((i < 0) or ((i >= len(_this))))) else _this[i])
                if ((char == "N") and multiposMap.h.get(i,None)):
                    line1.add("-1")
                else:
                    line1.add(SeqPhase1.code.h.get(char,None))
            _g6_head = varpos.h
            while (_g6_head is not None):
                val3 = _g6_head.item
                _g6_head = _g6_head.next
                i1 = val3
                _this1 = e2.getSeq()
                char1 = ("" if (((i1 < 0) or ((i1 >= len(_this1))))) else _this1[i1])
                if ((char1 == "N") and multiposMap.h.get(i1,None)):
                    line2.add("-1")
                else:
                    line2.add(SeqPhase1.code.h.get(char1,None))
            lines.add((HxOverrides.stringOrNull(line1.join(" ")) + " "))
            lines.add((HxOverrides.stringOrNull(line2.join(" ")) + " "))
        isOdd = False
        _g = 0
        _g1 = al2.getSequences()
        while (_g < len(_g1)):
            entry = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
            _g = (_g + 1)
            isOdd = (not isOdd)
            if isOdd:
                name = entry.getName()
                lines.add(HxString.substr(name,0,(len(name) - 1)))
            line = haxe_ds_List()
            _g5_head = varpos.h
            while (_g5_head is not None):
                val = _g5_head.item
                _g5_head = _g5_head.next
                i = val
                _this = entry.getSeq()
                char = ("" if (((i < 0) or ((i >= len(_this))))) else _this[i])
                if ((char == "N") and multiposMap.h.get(i,None)):
                    line.add("-1")
                else:
                    line.add(SeqPhase1.code.h.get(char,None))
            lines.add((HxOverrides.stringOrNull(line.join(" ")) + " "))
        isOdd = False
        _g = 0
        _g1 = al3.getSequences()
        while (_g < len(_g1)):
            entry = (_g1[_g] if _g >= 0 and _g < len(_g1) else None)
            _g = (_g + 1)
            isOdd = (not isOdd)
            if isOdd:
                name = entry.getName()
                lines.add(HxString.substr(name,0,(len(name) - 1)))
            line = haxe_ds_List()
            _g7_head = varpos.h
            while (_g7_head is not None):
                val = _g7_head.item
                _g7_head = _g7_head.next
                i = val
                _this = entry.getSeq()
                char = ("" if (((i < 0) or ((i >= len(_this))))) else _this[i])
                if ((char == "N") and multiposMap.h.get(i,None)):
                    line.add("-1")
                else:
                    line.add(SeqPhase1.code.h.get(char,None))
            lines.add((HxOverrides.stringOrNull(line.join(" ")) + " "))
        lines.add("")
        SeqPhase1Result.instance().setInpFile(lines.join("\n"))
        knownLines = haxe_ds_List()
        i = varpos.length
        result = haxe_ds_List()
        _g = 0
        _g1 = i
        while (_g < _g1):
            nnn = _g
            _g = (_g + 1)
            result.add("*")
        nStr = result.join("")
        i = varpos.length
        result = haxe_ds_List()
        _g = 0
        _g1 = i
        while (_g < _g1):
            nnn = _g
            _g = (_g + 1)
            result.add("0")
        oStr = result.join("")
        _g = 0
        _g1 = len(al1.getSequences())
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            knownLines.add(nStr)
        x = (len(al2.getSequences()) / 2)
        lll1 = None
        try:
            lll1 = int(x)
        except BaseException as _g:
            None
            lll1 = None
        _g = 0
        _g1 = lll1
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            knownLines.add(nStr)
        x = (len(al3.getSequences()) / 2)
        lll2 = None
        try:
            lll2 = int(x)
        except BaseException as _g:
            None
            lll2 = None
        _g = 0
        _g1 = lll2
        while (_g < _g1):
            i = _g
            _g = (_g + 1)
            knownLines.add(oStr)
        SeqPhase1Result.instance().setKnownFile(knownLines.join("\n"))
        if (len(al3.getSequences()) == 0):
            if (multipos.length == 0):
                SeqPhase1Result.instance().setSuggestedPhaseCommand("PHASE seqphase.inp seqphase.out")
            else:
                SeqPhase1Result.instance().setSuggestedPhaseCommand("PHASE -d1 seqphase.inp seqphase.out")
        elif (multipos.length == 0):
            SeqPhase1Result.instance().setSuggestedPhaseCommand("PHASE -kseqphase.known seqphase.inp seqphase.out")
        else:
            SeqPhase1Result.instance().setSuggestedPhaseCommand("PHASE -d1 -kseqphase.known seqphase.inp seqphase.out")
        return SeqPhase1Result.instance()

    @staticmethod
    def makeStr(c,i):
        result = haxe_ds_List()
        _g = 0
        _g1 = i
        while (_g < _g1):
            nnn = _g
            _g = (_g + 1)
            result.add(c)
        return result.join("")

    @staticmethod
    def main():
        Sys.stderr().writeString("SeqPHASE command-line version, Step 1: generating PHASE input files from FASTA alignments\n\n")
        Sys.stderr().writeString("Reference:\nFlot, J.-F. (2010) SeqPHASE: a web tool for interconverting PHASE input/output files and FASTA sequence alignments\nMolecular Ecology Ressources 10 (1): 162-166\n\n")
        Sys.stderr().writeString("Usage: perl seqphase1.pl -1 <first type of input file> -2 <second type of input file> -3 <third type of input file> -p <prefix>\nFirst type of input file = FASTA alignment of sequences from homozygous individuals and from heterozygotes to be phased (1 sequence per individual).\nSecond type of input file = FASTA alignment of alignment of fake haplotypes from heterozygotes to be phased (2 sequences per individual).\nThird type of input file (optional)= FASTA alignment of known haplotypes of previously phased heterozygotes (2 sequences per individual).\nPrefix for output files is optional, default prefix is 'phase'.\n\n")
        myArgs = Sys.args()
        file1 = None
        file2 = None
        file3 = None
        prefix = "phase"
        i = 0
        while (i < len(myArgs)):
            current = (myArgs[i] if i >= 0 and i < len(myArgs) else None)
            if (current == "-1"):
                i = (i + 1)
                myArgs1 = i
                file1 = sys_io_File.getContent((myArgs[myArgs1] if myArgs1 >= 0 and myArgs1 < len(myArgs) else None))
            elif (current == "-2"):
                i = (i + 1)
                myArgs2 = i
                file2 = sys_io_File.getContent((myArgs[myArgs2] if myArgs2 >= 0 and myArgs2 < len(myArgs) else None))
            elif (current == "-3"):
                i = (i + 1)
                myArgs3 = i
                file3 = sys_io_File.getContent((myArgs[myArgs3] if myArgs3 >= 0 and myArgs3 < len(myArgs) else None))
            elif (current == "-p"):
                i = (i + 1)
                prefix1 = i
                prefix = (myArgs[prefix1] if prefix1 >= 0 and prefix1 < len(myArgs) else None)
            else:
                Sys.stderr().writeString((("Error: Unknown commandline option " + ("null" if current is None else current)) + "\n"))
                Sys.exit(1)
            i = (i + 1)
        if ((file1 is None) and ((file2 is None))):
            Sys.stderr().writeString("There is no input file to be phased!\n")
            Sys.exit(1)
        result = SeqPhase1.doIt(file1,file2,file3)
        if result.hasAlign1Errors():
            Sys.stderr().writeString("The following errors occured while processing the -1 file:\n")
            _g = 0
            _g1 = result.getAlign1Errors()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasAlign1Warn():
            Sys.stderr().writeString("The following warnings occured while processing the -1 file:\n")
            _g = 0
            _g1 = result.getAlign1Warn()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasAlign2Errors():
            Sys.stderr().writeString("The following errors occured while processing the -2 file:\n")
            _g = 0
            _g1 = result.getAlign2Errors()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasAlign2Warn():
            Sys.stderr().writeString("The following warnings occured while processing the -2 file:\n")
            _g = 0
            _g1 = result.getAlign2Warn()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasAlign3Errors():
            Sys.stderr().writeString("The following errors occured while processing the -3 file:\n")
            _g = 0
            _g1 = result.getAlign3Errors()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasAlign3Warn():
            Sys.stderr().writeString("The following warnings occured while processing the -3 file:\n")
            _g = 0
            _g1 = result.getAlign3Warn()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasGeneralErrors():
            Sys.stderr().writeString("The following errors occured:\n")
            _g = 0
            _g1 = result.getGeneralErrors()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasGeneralWarn():
            Sys.stderr().writeString("The following warnings occured:\n")
            _g = 0
            _g1 = result.getGeneralWarn()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasNotes():
            _g = 0
            _g1 = result.getNotes()
            while (_g < len(_g1)):
                s = _g1[_g]
                _g = (_g + 1)
                Sys.stderr().writeString(s)
                Sys.stderr().writeString("\n")
        if result.hasInpFile():
            fout = sys_io_File.write((("null" if prefix is None else prefix) + ".inp"),False)
            fout.writeString(result.getInpFile())
            fout.close()
            Sys.stderr().writeString((("Since PHASE only accepts numbers and not letters for nucleotides at multistate characters, ? and N (missing information) have been replaced with ? or -1 (depending on whether the position displays two or more than two different nucleotides), - with 0, A with 1, C with 2, G with 3 and T with 4 in " + ("null" if prefix is None else prefix)) + ".inp (the main PHASE input file).\n"))
        if result.hasKnownFile():
            fout = sys_io_File.write((("null" if prefix is None else prefix) + ".known"),False)
            fout.writeString(result.getKnownFile())
            fout.close()
            Sys.stderr().writeString((("A known phase file (" + ("null" if prefix is None else prefix)) + ".known) has also been generated to tell PHASE which phases are already known and which ones are to infer.\n"))
        if result.hasConstFile():
            fout = sys_io_File.write((("null" if prefix is None else prefix) + ".const"),False)
            fout.writeString(result.getConstFile())
            fout.close()
            Sys.stderr().writeString((("In order to reduce PHASE running time, constant positions have been removed from the dataset. These constant positions have been stored in " + ("null" if prefix is None else prefix)) + ".const interspersed with periods (.) representing variable positions.\n"))
        if result.hasSuggestedCommand():
            Sys.stderr().writeString("Suggested command syntax (PHASE v2.1):\n\n")
            r = result.getSuggestedPhaseCommand()
            _this = r.split("seqphase")
            r = prefix.join([python_Boot.toString1(x1,'') for x1 in _this])
            Sys.stderr().writeString(r)
            Sys.stderr().writeString("\n\n")
            Sys.stderr().writeString((("(the result of the analysis will be stored in a series of files with names starting with '" + ("null" if prefix is None else prefix)) + ".out').\n"))
        Sys.stderr().writeString("\n")


class SeqPhase1Result:
    _hx_class_name = "SeqPhase1Result"
    __slots__ = ("errorsAlign1", "warningsAlign1", "errorsAlign2", "warningsAlign2", "errorsAlign3", "warningsAlign3", "errorsGeneral", "warningsGeneral", "notes", "suggestedPhaseCommand", "varPos", "nbVarPos", "inpFile", "knownFile", "constFile")
    _hx_fields = ["errorsAlign1", "warningsAlign1", "errorsAlign2", "warningsAlign2", "errorsAlign3", "warningsAlign3", "errorsGeneral", "warningsGeneral", "notes", "suggestedPhaseCommand", "varPos", "nbVarPos", "inpFile", "knownFile", "constFile"]
    _hx_methods = ["clear", "addErr", "addWrn", "addAlign1Error", "hasAlign1Errors", "getAlign1Errors", "addAlign1Warn", "hasAlign1Warn", "getAlign1Warn", "addAlign2Error", "hasAlign2Errors", "getAlign2Errors", "addAlign2Warn", "hasAlign2Warn", "getAlign2Warn", "addAlign3Error", "hasAlign3Errors", "getAlign3Errors", "addAlign3Warn", "hasAlign3Warn", "getAlign3Warn", "addGeneralError", "hasGeneralErrors", "getGeneralErrors", "addGeneralWarn", "hasGeneralWarn", "getGeneralWarn", "hasErrors", "addNote", "hasNotes", "getNotes", "setSuggestedPhaseCommand", "hasSuggestedCommand", "getSuggestedPhaseCommand", "setNrVarPos", "getNrVarPos", "setNrNbVarPos", "getNrNbVarPos", "setInpFile", "hasInpFile", "getInpFile", "setKnownFile", "hasKnownFile", "getKnownFile", "setConstFile", "hasConstFile", "getConstFile"]
    _hx_statics = ["inst", "instance"]

    def __init__(self):
        self.constFile = None
        self.knownFile = None
        self.inpFile = None
        self.nbVarPos = None
        self.varPos = None
        self.suggestedPhaseCommand = None
        self.notes = None
        self.warningsGeneral = None
        self.errorsGeneral = None
        self.warningsAlign3 = None
        self.errorsAlign3 = None
        self.warningsAlign2 = None
        self.errorsAlign2 = None
        self.warningsAlign1 = None
        self.errorsAlign1 = None
        self.clear()

    def clear(self):
        self.errorsAlign1 = haxe_ds_List()
        self.warningsAlign1 = haxe_ds_List()
        self.errorsAlign2 = haxe_ds_List()
        self.warningsAlign2 = haxe_ds_List()
        self.errorsAlign3 = haxe_ds_List()
        self.warningsAlign3 = haxe_ds_List()
        self.errorsGeneral = haxe_ds_List()
        self.warningsGeneral = haxe_ds_List()
        self.notes = haxe_ds_List()
        self.suggestedPhaseCommand = None
        self.varPos = None
        self.nbVarPos = None
        self.inpFile = None
        self.knownFile = None
        self.constFile = None

    def addErr(self,err,nr):
        if (nr == 1):
            self.addAlign1Error(err)
        elif (nr == 2):
            self.addAlign2Error(err)
        elif (nr == 3):
            self.addAlign3Error(err)
        else:
            raise haxe_Exception.thrown(("Illegal nr " + Std.string(nr)))

    def addWrn(self,wrn,nr):
        if (nr == 1):
            self.addAlign1Warn(wrn)
        elif (nr == 2):
            self.addAlign2Warn(wrn)
        elif (nr == 3):
            self.addAlign3Warn(wrn)
        else:
            raise haxe_Exception.thrown(("Illegal nr " + Std.string(nr)))

    def addAlign1Error(self,err):
        self.errorsAlign1.add(err)

    def hasAlign1Errors(self):
        return (self.errorsAlign1.length != 0)

    def getAlign1Errors(self):
        this1 = [None]*self.errorsAlign1.length
        result = this1
        i = 0
        _g_head = self.errorsAlign1.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def addAlign1Warn(self,wrn):
        self.warningsAlign1.add(wrn)

    def hasAlign1Warn(self):
        return (self.warningsAlign1.length != 0)

    def getAlign1Warn(self):
        this1 = [None]*self.warningsAlign1.length
        result = this1
        i = 0
        _g_head = self.warningsAlign1.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def addAlign2Error(self,err):
        self.errorsAlign2.add(err)

    def hasAlign2Errors(self):
        return (self.errorsAlign2.length != 0)

    def getAlign2Errors(self):
        this1 = [None]*self.errorsAlign2.length
        result = this1
        i = 0
        _g_head = self.errorsAlign2.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def addAlign2Warn(self,wrn):
        self.warningsAlign2.add(wrn)

    def hasAlign2Warn(self):
        return (self.warningsAlign2.length != 0)

    def getAlign2Warn(self):
        this1 = [None]*self.warningsAlign2.length
        result = this1
        i = 0
        _g_head = self.warningsAlign2.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def addAlign3Error(self,err):
        self.errorsAlign3.add(err)

    def hasAlign3Errors(self):
        return (self.errorsAlign3.length != 0)

    def getAlign3Errors(self):
        this1 = [None]*self.errorsAlign3.length
        result = this1
        i = 0
        _g_head = self.errorsAlign3.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def addAlign3Warn(self,wrn):
        self.warningsAlign3.add(wrn)

    def hasAlign3Warn(self):
        return (self.warningsAlign3.length != 0)

    def getAlign3Warn(self):
        this1 = [None]*self.warningsAlign3.length
        result = this1
        i = 0
        _g_head = self.warningsAlign3.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def addGeneralError(self,err):
        self.errorsGeneral.add(err)

    def hasGeneralErrors(self):
        return (self.errorsGeneral.length != 0)

    def getGeneralErrors(self):
        this1 = [None]*self.errorsGeneral.length
        result = this1
        i = 0
        _g_head = self.errorsGeneral.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def addGeneralWarn(self,wrn):
        self.warningsGeneral.add(wrn)

    def hasGeneralWarn(self):
        return (self.warningsGeneral.length != 0)

    def getGeneralWarn(self):
        this1 = [None]*self.warningsGeneral.length
        result = this1
        i = 0
        _g_head = self.warningsGeneral.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def hasErrors(self):
        if (not ((((self.errorsAlign1.length > 0) or ((self.errorsAlign2.length > 0))) or ((self.errorsAlign3.length > 0))))):
            return (self.errorsGeneral.length > 0)
        else:
            return True

    def addNote(self,note):
        self.notes.add(note)

    def hasNotes(self):
        return (self.notes.length != 0)

    def getNotes(self):
        this1 = [None]*self.notes.length
        result = this1
        i = 0
        _g_head = self.notes.h
        while (_g_head is not None):
            val = _g_head.item
            _g_head = _g_head.next
            item = val
            index = i
            i = (i + 1)
            result[index] = item
        return result

    def setSuggestedPhaseCommand(self,ph):
        self.suggestedPhaseCommand = ph

    def hasSuggestedCommand(self):
        return (self.suggestedPhaseCommand is not None)

    def getSuggestedPhaseCommand(self):
        return self.suggestedPhaseCommand

    def setNrVarPos(self,nr):
        self.varPos = nr

    def getNrVarPos(self):
        return self.varPos

    def setNrNbVarPos(self,nr):
        self.nbVarPos = nr

    def getNrNbVarPos(self):
        return self.nbVarPos

    def setInpFile(self,content):
        self.inpFile = content

    def hasInpFile(self):
        return (self.inpFile is not None)

    def getInpFile(self):
        return self.inpFile

    def setKnownFile(self,content):
        self.knownFile = content

    def hasKnownFile(self):
        return (self.knownFile is not None)

    def getKnownFile(self):
        return self.knownFile

    def setConstFile(self,content):
        self.constFile = content

    def hasConstFile(self):
        return (self.constFile is not None)

    def getConstFile(self):
        return self.constFile
    inst = None

    @staticmethod
    def instance():
        if (SeqPhase1Result.inst is None):
            SeqPhase1Result.inst = SeqPhase1Result()
        return SeqPhase1Result.inst



class Std:
    _hx_class_name = "Std"
    __slots__ = ()
    _hx_statics = ["isOfType", "string"]

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


class Float: pass


class Int: pass


class Bool: pass


class Dynamic: pass


class StringTools:
    _hx_class_name = "StringTools"
    __slots__ = ()
    _hx_statics = ["replace"]

    @staticmethod
    def replace(s,sub,by):
        _this = (list(s) if ((sub == "")) else s.split(sub))
        return by.join([python_Boot.toString1(x1,'') for x1 in _this])


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
    _hx_methods = ["add", "join"]

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
    _hx_methods = ["keys"]
    _hx_interfaces = [haxe_IMap]

    def __init__(self):
        self.h = dict()

    def keys(self):
        return python_HaxeIterator(iter(self.h.keys()))



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
    _hx_statics = ["eq", "stringOrNull", "modf", "mod"]

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
    def modf(a,b):
        if (b == 0.0):
            return float('nan')
        elif (a < 0):
            if (b < 0):
                return -(-a % (-b))
            else:
                return -(-a % b)
        elif (b < 0):
            return a % (-b)
        else:
            return a % b

    @staticmethod
    def mod(a,b):
        if (a < 0):
            if (b < 0):
                return -(-a % (-b))
            else:
                return -(-a % b)
        elif (b < 0):
            return a % (-b)
        else:
            return a % b


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

def _hx_init_FastaAlignmentParser_authorizedCharacters():
    def _hx_local_0():
        _g = haxe_ds_StringMap()
        _g.h["A"] = True
        _g.h["T"] = True
        _g.h["G"] = True
        _g.h["C"] = True
        _g.h["N"] = True
        _g.h["-"] = True
        _g.h["?"] = True
        _g.h["R"] = False
        _g.h["Y"] = False
        _g.h["M"] = False
        _g.h["K"] = False
        _g.h["W"] = False
        _g.h["S"] = False
        return _g
    return _hx_local_0()
FastaAlignmentParser.authorizedCharacters = _hx_init_FastaAlignmentParser_authorizedCharacters()
def _hx_init_SeqPhase1_map1():
    def _hx_local_0():
        _g = haxe_ds_StringMap()
        _g.h["W"] = "A"
        _g.h["S"] = "C"
        _g.h["K"] = "T"
        _g.h["M"] = "A"
        _g.h["Y"] = "C"
        _g.h["R"] = "A"
        return _g
    return _hx_local_0()
SeqPhase1.map1 = _hx_init_SeqPhase1_map1()
def _hx_init_SeqPhase1_map2():
    def _hx_local_0():
        _g = haxe_ds_StringMap()
        _g.h["W"] = "T"
        _g.h["S"] = "G"
        _g.h["K"] = "G"
        _g.h["M"] = "C"
        _g.h["Y"] = "T"
        _g.h["R"] = "G"
        return _g
    return _hx_local_0()
SeqPhase1.map2 = _hx_init_SeqPhase1_map2()
def _hx_init_SeqPhase1_code():
    def _hx_local_0():
        _g = haxe_ds_StringMap()
        _g.h["A"] = "1"
        _g.h["C"] = "2"
        _g.h["G"] = "3"
        _g.h["T"] = "4"
        _g.h["?"] = "?"
        _g.h["N"] = "?"
        _g.h["-"] = "0"
        return _g
    return _hx_local_0()
SeqPhase1.code = _hx_init_SeqPhase1_code()
python_Boot.keywords = set(["and", "del", "from", "not", "with", "as", "elif", "global", "or", "yield", "assert", "else", "if", "pass", "None", "break", "except", "import", "raise", "True", "class", "exec", "in", "return", "False", "continue", "finally", "is", "try", "def", "for", "lambda", "while"])
python_Boot.prefixLength = len("_hx_")

SeqPhase1.main()
