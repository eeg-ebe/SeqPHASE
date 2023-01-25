function $extend(from, fields) {
	var proto = Object.create(from);
	for (var name in fields) proto[name] = fields[name];
	if( fields.toString !== Object.prototype.toString ) proto.toString = fields.toString;
	return proto;
}
var Entry = function(line,name,seq) {
	if(line == null) {
		line = -1;
	}
	this.line = line;
	this.name = name;
	this.seq = seq;
};
Entry.__name__ = true;
Entry.prototype = {
	addToSeq: function(s) {
		if(this.seq == null || this.seq == "") {
			this.seq = s;
		} else {
			this.seq += s;
		}
	}
	,getName: function() {
		return this.name;
	}
	,getSeq: function() {
		if(this.seq == null) {
			return "";
		}
		return this.seq;
	}
	,getLineNo: function() {
		return this.line;
	}
};
var FastaAlignmentParser = function(fileContent,allChecks,allSort,fileNr) {
	this.fastaContent = [];
	this.seqLength = -1;
	if(fileContent == null) {
		return;
	}
	var end = fileContent.length;
	while(true) {
		var tmp;
		if(end > 0) {
			var cCode = HxOverrides.cca(fileContent,end - 1);
			var result = false;
			var _g = 0;
			var _g1 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
			while(_g < _g1.length) {
				var ele = _g1[_g];
				++_g;
				if(ele == cCode) {
					result = true;
					break;
				}
			}
			tmp = result;
		} else {
			tmp = false;
		}
		if(!tmp) {
			break;
		}
		--end;
	}
	var s = fileContent.substring(0,end);
	var begin = 0;
	var sLen = s.length;
	while(true) {
		var tmp;
		if(begin < sLen) {
			var cCode = HxOverrides.cca(s,begin);
			var result = false;
			var _g = 0;
			var _g1 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
			while(_g < _g1.length) {
				var ele = _g1[_g];
				++_g;
				if(ele == cCode) {
					result = true;
					break;
				}
			}
			tmp = result;
		} else {
			tmp = false;
		}
		if(!tmp) {
			break;
		}
		++begin;
	}
	if(HxOverrides.substr(s,begin,null) == "") {
		SeqPhase1Result.instance().addErr("Empty file!",fileNr);
		return;
	}
	if(HxOverrides.substr(fileContent,0,">".length) != ">" && HxOverrides.substr(fileContent,0,";".length) != ";") {
		SeqPhase1Result.instance().addErr("File does not seem to be a fasta file!",fileNr);
		return;
	}
	var lines = fileContent.split("\n");
	if(lines.length == 0) {
		SeqPhase1Result.instance().addErr("Not a fasta but an empty file!",fileNr);
		return;
	} else if(lines.length == 1) {
		SeqPhase1Result.instance().addErr("Only 1 line detected! Please check data format (opening the alignment in MEGA (http://www.megasoftware.net/) and exporting it as FASTA again may solve the problem; alternatively, there may be a problem with end-of-line characters - see http://en.wikipedia.org/wiki/Newline for details!",fileNr);
		return;
	}
	var entryMap_h = Object.create(null);
	var current = null;
	var lineNo = 0;
	var underscoreWarningOutputted = false;
	var _g = 0;
	while(_g < lines.length) {
		var line = lines[_g];
		++_g;
		++lineNo;
		var end = line.length;
		while(true) {
			var line1;
			if(end > 0) {
				var cCode = HxOverrides.cca(line,end - 1);
				var result = false;
				var _g1 = 0;
				var _g2 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
				while(_g1 < _g2.length) {
					var ele = _g2[_g1];
					++_g1;
					if(ele == cCode) {
						result = true;
						break;
					}
				}
				line1 = result;
			} else {
				line1 = false;
			}
			if(!line1) {
				break;
			}
			--end;
		}
		var s = line.substring(0,end);
		var begin = 0;
		var sLen = s.length;
		while(true) {
			var line2;
			if(begin < sLen) {
				var cCode1 = HxOverrides.cca(s,begin);
				var result1 = false;
				var _g3 = 0;
				var _g4 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
				while(_g3 < _g4.length) {
					var ele1 = _g4[_g3];
					++_g3;
					if(ele1 == cCode1) {
						result1 = true;
						break;
					}
				}
				line2 = result1;
			} else {
				line2 = false;
			}
			if(!line2) {
				break;
			}
			++begin;
		}
		line = HxOverrides.substr(s,begin,null);
		if(HxOverrides.substr(line,0,";".length) == ";") {
			continue;
		} else if(HxOverrides.substr(line,0,">".length) == ">") {
			var s1 = HxOverrides.substr(line,1,null);
			var end1 = s1.length;
			while(true) {
				var tmp;
				if(end1 > 0) {
					var cCode2 = HxOverrides.cca(s1,end1 - 1);
					var result2 = false;
					var _g5 = 0;
					var _g6 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
					while(_g5 < _g6.length) {
						var ele2 = _g6[_g5];
						++_g5;
						if(ele2 == cCode2) {
							result2 = true;
							break;
						}
					}
					tmp = result2;
				} else {
					tmp = false;
				}
				if(!tmp) {
					break;
				}
				--end1;
			}
			var s2 = s1.substring(0,end1);
			var begin1 = 0;
			var sLen1 = s2.length;
			while(true) {
				var tmp1;
				if(begin1 < sLen1) {
					var cCode3 = HxOverrides.cca(s2,begin1);
					var result3 = false;
					var _g7 = 0;
					var _g8 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
					while(_g7 < _g8.length) {
						var ele3 = _g8[_g7];
						++_g7;
						if(ele3 == cCode3) {
							result3 = true;
							break;
						}
					}
					tmp1 = result3;
				} else {
					tmp1 = false;
				}
				if(!tmp1) {
					break;
				}
				++begin1;
			}
			var indName = HxOverrides.substr(s2,begin1,null);
			var indNameCor = StringTools.replace(indName," ","_");
			if(indName != indNameCor) {
				if(!underscoreWarningOutputted) {
					SeqPhase1Result.instance().addWrn("Warning: PHASE does not accept spaces in individual names. These spaces got replaced by underscore characters.",fileNr);
					underscoreWarningOutputted = true;
				}
				indName = indNameCor;
			}
			if(Object.prototype.hasOwnProperty.call(entryMap_h,indName)) {
				SeqPhase1Result.instance().addErr("Repeat of name " + indName + " encountered in alignment (line " + lineNo + ", line " + entryMap_h[indName].getLineNo() + ")",fileNr);
			}
			current = new Entry(lineNo,indName);
			entryMap_h[indName] = current;
			if(indName.length == 0) {
				SeqPhase1Result.instance().addErr("Missing sequence name, line " + lineNo,fileNr);
			}
		} else {
			var _g9 = 0;
			var _g10 = line.length;
			while(_g9 < _g10) {
				var i = _g9++;
				var char = line.charAt(i).toUpperCase();
				if(!Object.prototype.hasOwnProperty.call(FastaAlignmentParser.authorizedCharacters.h,char)) {
					SeqPhase1Result.instance().addErr("Unknown character state " + char + " in " + current.getName() + ", line " + lineNo,fileNr);
				} else if(allChecks && FastaAlignmentParser.authorizedCharacters.h[char] == false) {
					SeqPhase1Result.instance().addErr("Unallowed state " + char + " in " + current.getName() + ", line " + lineNo,fileNr);
				}
			}
			line = line.split("?").join("N");
			current.addToSeq(line.toUpperCase());
		}
	}
	if(current == null) {
		SeqPhase1Result.instance().addErr("Corrupted Fasta File",fileNr);
		return;
	}
	if(current.getSeq() == "") {
		SeqPhase1Result.instance().addErr("Empty sequence " + current.getName() + ", line " + current.getLineNo(),fileNr);
		return;
	}
	this.seqLength = current.getSeq().length;
	var h = entryMap_h;
	var _g1_h = h;
	var _g1_keys = Object.keys(h);
	var _g1_length = _g1_keys.length;
	var _g1_current = 0;
	while(_g1_current < _g1_length) {
		var key = _g1_keys[_g1_current++];
		var val = entryMap_h[key];
		this.fastaContent.push(val);
		if(val.getSeq().length != current.getSeq().length) {
			SeqPhase1Result.instance().addErr("Not all sequences in this file have equal lengths. E.g. sequence " + val.getName() + " (line " + val.getLineNo() + ") is of length " + val.getSeq().length + " while sequence " + current.getName() + " (line " + current.getLineNo() + ") is of length " + current.getSeq().length,fileNr);
			return;
		}
		if(HxOverrides.substr(val.getSeq(),0,"-".length) == "-") {
			SeqPhase1Result.instance().addWrn("Sequence " + val.getName() + " (line " + val.getLineNo() + ") starts with '-'. Is it a real indel or did you mean 'N' or '?' (missing data)?",fileNr);
		}
		if(val.getSeq().charAt(val.getSeq().length - 1) == "-") {
			SeqPhase1Result.instance().addWrn("Sequence " + val.getName() + " (line " + val.getLineNo() + ") ends with '-'. Is it a real indel or did you mean 'N' or '?' (missing data)?",fileNr);
		}
	}
	this.fastaContent.sort(function(e1,e2) {
		if(allSort) {
			var nameE1 = e1.getName();
			var nameE2 = e1.getName();
			var aNameE1 = nameE1.substring(0,nameE1.length - 1);
			var aNameE2 = nameE2.substring(0,nameE2.length - 1);
			if(aNameE1 < aNameE2) {
				return -1;
			}
			if(aNameE1 > aNameE2) {
				return 1;
			}
			var lNameE1 = nameE1.charAt(nameE1.length - 1);
			var lNameE2 = nameE2.charAt(nameE2.length - 1);
			if(lNameE1 < lNameE2) {
				return -1;
			}
			if(lNameE1 > lNameE2) {
				return 1;
			}
			return 0;
		} else {
			if(e1.getName() < e2.getName()) {
				return -1;
			}
			if(e1.getName() > e2.getName()) {
				return 1;
			}
			return 0;
		}
	});
	if(allChecks) {
		if(this.fastaContent.length % 2 != 0) {
			SeqPhase1Result.instance().addErr("Uneven number of sequences in alignment: please check data.",fileNr);
		}
		var lastName = null;
		var _g = 0;
		var _g1 = this.fastaContent;
		while(_g < _g1.length) {
			var entry = _g1[_g];
			++_g;
			if(lastName == null) {
				lastName = entry.getName();
				lastName = lastName.substring(0,lastName.length - 1);
			} else {
				var curName = entry.getName();
				curName = curName.substring(0,curName.length - 1);
				if(lastName == curName) {
					lastName = null;
				}
			}
		}
		var map_h = Object.create(null);
		var _g = 0;
		var _g1 = this.fastaContent;
		while(_g < _g1.length) {
			var entry = _g1[_g];
			++_g;
			var name = entry.getName();
			name = name.substring(0,name.length - 1);
			if(Object.prototype.hasOwnProperty.call(map_h,name)) {
				var i = map_h[name];
				map_h[name] = i + 1;
			} else {
				map_h[name] = 1;
			}
		}
		var h = map_h;
		var _g6_h = h;
		var _g6_keys = Object.keys(h);
		var _g6_length = _g6_keys.length;
		var _g6_current = 0;
		while(_g6_current < _g6_length) {
			var name = _g6_keys[_g6_current++];
			var i = map_h[name];
			if(i != 2) {
				var lst = new haxe_ds_List();
				var _g = 0;
				var _g1 = this.fastaContent;
				while(_g < _g1.length) {
					var entry = _g1[_g];
					++_g;
					var nameEntry = entry.getName();
					var subnameEntry = nameEntry.substring(0,nameEntry.length - 1);
					if(subnameEntry == name) {
						lst.add(nameEntry);
					}
				}
				if(i == 1) {
					SeqPhase1Result.instance().addErr("Found only one haplotype sequence (" + lst.join(",") + ") for individual '" + name + "'!",fileNr);
				} else {
					SeqPhase1Result.instance().addErr("Found " + i + " haplotype sequences (" + lst.join(",") + ") for individual '" + name + "'!",fileNr);
				}
			}
		}
	}
};
FastaAlignmentParser.__name__ = true;
FastaAlignmentParser.startsWith = function(t,s) {
	return HxOverrides.substr(t,0,s.length) == s;
};
FastaAlignmentParser.isWhitespace = function(s,pos) {
	var cCode = HxOverrides.cca(s,pos);
	var result = false;
	var _g = 0;
	var _g1 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
	while(_g < _g1.length) {
		var ele = _g1[_g];
		++_g;
		if(ele == cCode) {
			result = true;
			break;
		}
	}
	return result;
};
FastaAlignmentParser.stripStringBegin = function(s) {
	var begin = 0;
	var sLen = s.length;
	while(true) {
		var tmp;
		if(begin < sLen) {
			var cCode = HxOverrides.cca(s,begin);
			var result = false;
			var _g = 0;
			var _g1 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
			while(_g < _g1.length) {
				var ele = _g1[_g];
				++_g;
				if(ele == cCode) {
					result = true;
					break;
				}
			}
			tmp = result;
		} else {
			tmp = false;
		}
		if(!tmp) {
			break;
		}
		++begin;
	}
	return HxOverrides.substr(s,begin,null);
};
FastaAlignmentParser.stripStringEnd = function(s) {
	var end = s.length;
	while(true) {
		var tmp;
		if(end > 0) {
			var cCode = HxOverrides.cca(s,end - 1);
			var result = false;
			var _g = 0;
			var _g1 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
			while(_g < _g1.length) {
				var ele = _g1[_g];
				++_g;
				if(ele == cCode) {
					result = true;
					break;
				}
			}
			tmp = result;
		} else {
			tmp = false;
		}
		if(!tmp) {
			break;
		}
		--end;
	}
	return s.substring(0,end);
};
FastaAlignmentParser.stripString = function(s) {
	var end = s.length;
	while(true) {
		var tmp;
		if(end > 0) {
			var cCode = HxOverrides.cca(s,end - 1);
			var result = false;
			var _g = 0;
			var _g1 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
			while(_g < _g1.length) {
				var ele = _g1[_g];
				++_g;
				if(ele == cCode) {
					result = true;
					break;
				}
			}
			tmp = result;
		} else {
			tmp = false;
		}
		if(!tmp) {
			break;
		}
		--end;
	}
	var s1 = s.substring(0,end);
	var begin = 0;
	var sLen = s1.length;
	while(true) {
		var tmp;
		if(begin < sLen) {
			var cCode = HxOverrides.cca(s1,begin);
			var result = false;
			var _g = 0;
			var _g1 = [9,10,11,12,13,32,133,160,5760,8192,8192,8193,8194,8195,8196,8197,8198,8199,8200,8201,8202,8232,8233,8239,8287,12288,6158,8203,8204,8205,8288,65279];
			while(_g < _g1.length) {
				var ele = _g1[_g];
				++_g;
				if(ele == cCode) {
					result = true;
					break;
				}
			}
			tmp = result;
		} else {
			tmp = false;
		}
		if(!tmp) {
			break;
		}
		++begin;
	}
	return HxOverrides.substr(s1,begin,null);
};
FastaAlignmentParser.prototype = {
	getSeqLength: function() {
		return this.seqLength;
	}
	,getSequences: function() {
		return this.fastaContent;
	}
};
var HxOverrides = function() { };
HxOverrides.__name__ = true;
HxOverrides.cca = function(s,index) {
	var x = s.charCodeAt(index);
	if(x != x) {
		return undefined;
	}
	return x;
};
HxOverrides.substr = function(s,pos,len) {
	if(len == null) {
		len = s.length;
	} else if(len < 0) {
		if(pos == 0) {
			len = s.length + len;
		} else {
			return "";
		}
	}
	return s.substr(pos,len);
};
HxOverrides.now = function() {
	return Date.now();
};
Math.__name__ = true;
var SeqPhase1 = function() { };
SeqPhase1.__name__ = true;
SeqPhase1.doIt = function(align1,align2,align3) {
	SeqPhase1Result.instance().clear();
	var al1 = new FastaAlignmentParser(align1,false,false,1);
	var al2 = new FastaAlignmentParser(align2,true,true,2);
	var al3 = new FastaAlignmentParser(align3,true,false,3);
	var expectedLength = al1.getSeqLength() > al2.getSeqLength() ? al1.getSeqLength() : al2.getSeqLength();
	if(expectedLength <= al3.getSeqLength()) {
		expectedLength = al3.getSeqLength();
	}
	var diffLength = false;
	if(al1.getSeqLength() != -1 && al1.getSeqLength() != expectedLength) {
		diffLength = true;
	}
	if(al2.getSeqLength() != -1 && al2.getSeqLength() != expectedLength) {
		diffLength = true;
	}
	if(al3.getSeqLength() != -1 && al3.getSeqLength() != expectedLength) {
		diffLength = true;
	}
	if(diffLength) {
		SeqPhase1Result.instance().addGeneralError("Not all input sequences have equal lengths, please check whether this is expected.");
	} else if(expectedLength == -1 || expectedLength == 0) {
		SeqPhase1Result.instance().addGeneralError("It seems that all given sequences are empty ...");
	}
	if(SeqPhase1Result.instance().hasErrors()) {
		return SeqPhase1Result.instance();
	}
	var al1a = new haxe_ds_List();
	var al1b = new haxe_ds_List();
	var _g = 0;
	var _g1 = al1.getSequences();
	while(_g < _g1.length) {
		var entry = _g1[_g];
		++_g;
		var seq1a = new haxe_ds_List();
		var seq1b = new haxe_ds_List();
		var _g2 = 0;
		var _g3 = entry.getSeq().length;
		while(_g2 < _g3) {
			var i = _g2++;
			var c = entry.getSeq().charAt(i);
			if(Object.prototype.hasOwnProperty.call(SeqPhase1.map1.h,c)) {
				seq1a.add(SeqPhase1.map1.h[c]);
				seq1b.add(SeqPhase1.map2.h[c]);
			} else {
				seq1a.add(c);
				seq1b.add(c);
			}
		}
		al1a.add(new Entry(entry.getLineNo(),entry.getName(),seq1a.join("")));
		al1b.add(new Entry(entry.getLineNo(),entry.getName(),seq1b.join("")));
	}
	var varpos = new haxe_ds_List();
	var multipos = new haxe_ds_List();
	var multiposMap_h = { };
	var constFileContent = new haxe_ds_List();
	var _g = 0;
	var _g1 = expectedLength;
	while(_g < _g1) {
		var i = _g++;
		var m_h = Object.create(null);
		var _g2_head = al1a.h;
		while(_g2_head != null) {
			var val = _g2_head.item;
			_g2_head = _g2_head.next;
			var entry = val;
			m_h[entry.getSeq().charAt(i)] = false;
		}
		var _g3_head = al1b.h;
		while(_g3_head != null) {
			var val1 = _g3_head.item;
			_g3_head = _g3_head.next;
			var entry1 = val1;
			m_h[entry1.getSeq().charAt(i)] = false;
		}
		var _g2 = 0;
		var _g3 = al2.getSequences();
		while(_g2 < _g3.length) {
			var entry2 = _g3[_g2];
			++_g2;
			m_h[entry2.getSeq().charAt(i)] = false;
		}
		var _g4 = 0;
		var _g5 = al3.getSequences();
		while(_g4 < _g5.length) {
			var entry3 = _g5[_g4];
			++_g4;
			m_h[entry3.getSeq().charAt(i)] = false;
		}
		var mapLen = 0;
		var mapLenWithoutNs = 0;
		var lastKey = null;
		var h = m_h;
		var _g8_h = h;
		var _g8_keys = Object.keys(h);
		var _g8_length = _g8_keys.length;
		var _g8_current = 0;
		while(_g8_current < _g8_length) {
			var key = _g8_keys[_g8_current++];
			lastKey = key;
			++mapLen;
			if(key != "N") {
				++mapLenWithoutNs;
			}
		}
		if(mapLen == 0) {
			SeqPhase1Result.instance().addGeneralError("Bug detected: There seems to be a bug in this implementation of SeqPHASE. Please contact the author so that they can fix this bug.");
			constFileContent.add("X");
		} else if(mapLen == 1) {
			if(lastKey == "N") {
				SeqPhase1Result.instance().addGeneralWarn("Found only N/? s at position " + (i + 1) + ".");
			} else if(lastKey == "-") {
				SeqPhase1Result.instance().addGeneralWarn("Found only -'s at position " + (i + 1) + ". This may indicate an alignment problem. Consider to recheck your input data!");
			}
			constFileContent.add(lastKey);
		} else {
			constFileContent.add(".");
			varpos.add(i);
			if(mapLenWithoutNs > 2) {
				multipos.add(i);
				multiposMap_h[i] = true;
			} else {
				multiposMap_h[i] = false;
			}
		}
	}
	if(varpos.length == 0) {
		SeqPhase1Result.instance().addGeneralError("Not a single variable position detected in dataset! Please check data.");
	} else {
		SeqPhase1Result.instance().addNote("There are " + varpos.length + " variable positions in your dataset, including " + multipos.length + " position(s) with more than two different states.");
	}
	SeqPhase1Result.instance().setConstFile(constFileContent.join(""));
	var lines = new haxe_ds_List();
	lines.add("" + (al1a.length + al2.getSequences().length / 2 + al3.getSequences().length / 2 | 0));
	lines.add("" + varpos.length);
	var l1 = new haxe_ds_List();
	var l2 = new haxe_ds_List();
	var _g4_head = varpos.h;
	while(_g4_head != null) {
		var val = _g4_head.item;
		_g4_head = _g4_head.next;
		var pos = val;
		l1.add("" + (pos + 1));
		if(multiposMap_h[pos]) {
			l2.add("M");
		} else {
			l2.add("S");
		}
	}
	lines.add("P " + l1.join(" "));
	lines.add(l2.join(" ") + " ");
	var it1_head = al1a.h;
	var it2_head = al1b.h;
	while(it1_head != null) {
		var val = it1_head.item;
		it1_head = it1_head.next;
		var e1 = val;
		var val1 = it2_head.item;
		it2_head = it2_head.next;
		var e2 = val1;
		lines.add(e1.getName());
		var line1 = new haxe_ds_List();
		var line2 = new haxe_ds_List();
		var _g5_head = varpos.h;
		while(_g5_head != null) {
			var val2 = _g5_head.item;
			_g5_head = _g5_head.next;
			var i = val2;
			var char = e1.getSeq().charAt(i);
			if(char == "N" && multiposMap_h[i]) {
				line1.add("-1");
			} else {
				line1.add(SeqPhase1.code.h[char]);
			}
		}
		var _g6_head = varpos.h;
		while(_g6_head != null) {
			var val3 = _g6_head.item;
			_g6_head = _g6_head.next;
			var i1 = val3;
			var char1 = e2.getSeq().charAt(i1);
			if(char1 == "N" && multiposMap_h[i1]) {
				line2.add("-1");
			} else {
				line2.add(SeqPhase1.code.h[char1]);
			}
		}
		lines.add(line1.join(" ") + " ");
		lines.add(line2.join(" ") + " ");
	}
	var isOdd = false;
	var _g = 0;
	var _g1 = al2.getSequences();
	while(_g < _g1.length) {
		var entry = _g1[_g];
		++_g;
		isOdd = !isOdd;
		if(isOdd) {
			var name = entry.getName();
			lines.add(HxOverrides.substr(name,0,name.length - 1));
		}
		var line = new haxe_ds_List();
		var _g5_head = varpos.h;
		while(_g5_head != null) {
			var val = _g5_head.item;
			_g5_head = _g5_head.next;
			var i = val;
			var char = entry.getSeq().charAt(i);
			if(char == "N" && multiposMap_h[i]) {
				line.add("-1");
			} else {
				line.add(SeqPhase1.code.h[char]);
			}
		}
		lines.add(line.join(" ") + " ");
	}
	isOdd = false;
	var _g = 0;
	var _g1 = al3.getSequences();
	while(_g < _g1.length) {
		var entry = _g1[_g];
		++_g;
		isOdd = !isOdd;
		if(isOdd) {
			var name = entry.getName();
			lines.add(HxOverrides.substr(name,0,name.length - 1));
		}
		var line = new haxe_ds_List();
		var _g7_head = varpos.h;
		while(_g7_head != null) {
			var val = _g7_head.item;
			_g7_head = _g7_head.next;
			var i = val;
			var char = entry.getSeq().charAt(i);
			if(char == "N" && multiposMap_h[i]) {
				line.add("-1");
			} else {
				line.add(SeqPhase1.code.h[char]);
			}
		}
		lines.add(line.join(" ") + " ");
	}
	lines.add("");
	SeqPhase1Result.instance().setInpFile(lines.join("\n"));
	var knownLines = new haxe_ds_List();
	var i = varpos.length;
	var result = new haxe_ds_List();
	var _g = 0;
	var _g1 = i;
	while(_g < _g1) {
		var nnn = _g++;
		result.add("*");
	}
	var nStr = result.join("");
	var i = varpos.length;
	var result = new haxe_ds_List();
	var _g = 0;
	var _g1 = i;
	while(_g < _g1) {
		var nnn = _g++;
		result.add("0");
	}
	var oStr = result.join("");
	var _g = 0;
	var _g1 = al1.getSequences().length;
	while(_g < _g1) {
		var i = _g++;
		knownLines.add(nStr);
	}
	var lll1 = al2.getSequences().length / 2 | 0;
	var _g = 0;
	var _g1 = lll1;
	while(_g < _g1) {
		var i = _g++;
		knownLines.add(nStr);
	}
	var lll2 = al3.getSequences().length / 2 | 0;
	var _g = 0;
	var _g1 = lll2;
	while(_g < _g1) {
		var i = _g++;
		knownLines.add(oStr);
	}
	SeqPhase1Result.instance().setKnownFile(knownLines.join("\n"));
	if(al3.getSequences().length == 0) {
		if(multipos.length == 0) {
			SeqPhase1Result.instance().setSuggestedPhaseCommand("PHASE seqphase.inp seqphase.out");
		} else {
			SeqPhase1Result.instance().setSuggestedPhaseCommand("PHASE -d1 seqphase.inp seqphase.out");
		}
	} else if(multipos.length == 0) {
		SeqPhase1Result.instance().setSuggestedPhaseCommand("PHASE -kseqphase.known seqphase.inp seqphase.out");
	} else {
		SeqPhase1Result.instance().setSuggestedPhaseCommand("PHASE -d1 -kseqphase.known seqphase.inp seqphase.out");
	}
	return SeqPhase1Result.instance();
};
SeqPhase1.makeStr = function(c,i) {
	var result = new haxe_ds_List();
	var _g = 0;
	var _g1 = i;
	while(_g < _g1) {
		var nnn = _g++;
		result.add(c);
	}
	return result.join("");
};
SeqPhase1.main = function() {
};
var SeqPhase1Result = function() {
	this.clear();
};
SeqPhase1Result.__name__ = true;
SeqPhase1Result.instance = function() {
	if(SeqPhase1Result.inst == null) {
		SeqPhase1Result.inst = new SeqPhase1Result();
	}
	return SeqPhase1Result.inst;
};
SeqPhase1Result.prototype = {
	clear: function() {
		this.errorsAlign1 = new haxe_ds_List();
		this.warningsAlign1 = new haxe_ds_List();
		this.errorsAlign2 = new haxe_ds_List();
		this.warningsAlign2 = new haxe_ds_List();
		this.errorsAlign3 = new haxe_ds_List();
		this.warningsAlign3 = new haxe_ds_List();
		this.errorsGeneral = new haxe_ds_List();
		this.warningsGeneral = new haxe_ds_List();
		this.notes = new haxe_ds_List();
		this.suggestedPhaseCommand = null;
		this.varPos = null;
		this.nbVarPos = null;
		this.inpFile = null;
		this.knownFile = null;
		this.constFile = null;
	}
	,addErr: function(err,nr) {
		if(nr == 1) {
			this.addAlign1Error(err);
		} else if(nr == 2) {
			this.addAlign2Error(err);
		} else if(nr == 3) {
			this.addAlign3Error(err);
		} else {
			throw haxe_Exception.thrown("Illegal nr " + nr);
		}
	}
	,addWrn: function(wrn,nr) {
		if(nr == 1) {
			this.addAlign1Warn(wrn);
		} else if(nr == 2) {
			this.addAlign2Warn(wrn);
		} else if(nr == 3) {
			this.addAlign3Warn(wrn);
		} else {
			throw haxe_Exception.thrown("Illegal nr " + nr);
		}
	}
	,addAlign1Error: function(err) {
		this.errorsAlign1.add(err);
	}
	,hasAlign1Errors: function() {
		return this.errorsAlign1.length != 0;
	}
	,getAlign1Errors: function() {
		var this1 = new Array(this.errorsAlign1.length);
		var result = this1;
		var i = 0;
		var _g_head = this.errorsAlign1.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,addAlign1Warn: function(wrn) {
		this.warningsAlign1.add(wrn);
	}
	,hasAlign1Warn: function() {
		return this.warningsAlign1.length != 0;
	}
	,getAlign1Warn: function() {
		var this1 = new Array(this.warningsAlign1.length);
		var result = this1;
		var i = 0;
		var _g_head = this.warningsAlign1.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,addAlign2Error: function(err) {
		this.errorsAlign2.add(err);
	}
	,hasAlign2Errors: function() {
		return this.errorsAlign2.length != 0;
	}
	,getAlign2Errors: function() {
		var this1 = new Array(this.errorsAlign2.length);
		var result = this1;
		var i = 0;
		var _g_head = this.errorsAlign2.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,addAlign2Warn: function(wrn) {
		this.warningsAlign2.add(wrn);
	}
	,hasAlign2Warn: function() {
		return this.warningsAlign2.length != 0;
	}
	,getAlign2Warn: function() {
		var this1 = new Array(this.warningsAlign2.length);
		var result = this1;
		var i = 0;
		var _g_head = this.warningsAlign2.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,addAlign3Error: function(err) {
		this.errorsAlign3.add(err);
	}
	,hasAlign3Errors: function() {
		return this.errorsAlign3.length != 0;
	}
	,getAlign3Errors: function() {
		var this1 = new Array(this.errorsAlign3.length);
		var result = this1;
		var i = 0;
		var _g_head = this.errorsAlign3.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,addAlign3Warn: function(wrn) {
		this.warningsAlign3.add(wrn);
	}
	,hasAlign3Warn: function() {
		return this.warningsAlign3.length != 0;
	}
	,getAlign3Warn: function() {
		var this1 = new Array(this.warningsAlign3.length);
		var result = this1;
		var i = 0;
		var _g_head = this.warningsAlign3.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,addGeneralError: function(err) {
		this.errorsGeneral.add(err);
	}
	,hasGeneralErrors: function() {
		return this.errorsGeneral.length != 0;
	}
	,getGeneralErrors: function() {
		var this1 = new Array(this.errorsGeneral.length);
		var result = this1;
		var i = 0;
		var _g_head = this.errorsGeneral.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,addGeneralWarn: function(wrn) {
		this.warningsGeneral.add(wrn);
	}
	,hasGeneralWarn: function() {
		return this.warningsGeneral.length != 0;
	}
	,getGeneralWarn: function() {
		var this1 = new Array(this.warningsGeneral.length);
		var result = this1;
		var i = 0;
		var _g_head = this.warningsGeneral.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,hasErrors: function() {
		if(!(this.errorsAlign1.length > 0 || this.errorsAlign2.length > 0 || this.errorsAlign3.length > 0)) {
			return this.errorsGeneral.length > 0;
		} else {
			return true;
		}
	}
	,addNote: function(note) {
		this.notes.add(note);
	}
	,hasNotes: function() {
		return this.notes.length != 0;
	}
	,getNotes: function() {
		var this1 = new Array(this.notes.length);
		var result = this1;
		var i = 0;
		var _g_head = this.notes.h;
		while(_g_head != null) {
			var val = _g_head.item;
			_g_head = _g_head.next;
			var item = val;
			result[i++] = item;
		}
		return result;
	}
	,setSuggestedPhaseCommand: function(ph) {
		this.suggestedPhaseCommand = ph;
	}
	,hasSuggestedCommand: function() {
		return this.suggestedPhaseCommand != null;
	}
	,getSuggestedPhaseCommand: function() {
		return this.suggestedPhaseCommand;
	}
	,setNrVarPos: function(nr) {
		this.varPos = nr;
	}
	,getNrVarPos: function() {
		return this.varPos;
	}
	,setNrNbVarPos: function(nr) {
		this.nbVarPos = nr;
	}
	,getNrNbVarPos: function() {
		return this.nbVarPos;
	}
	,setInpFile: function(content) {
		this.inpFile = content;
	}
	,hasInpFile: function() {
		return this.inpFile != null;
	}
	,getInpFile: function() {
		return this.inpFile;
	}
	,setKnownFile: function(content) {
		this.knownFile = content;
	}
	,hasKnownFile: function() {
		return this.knownFile != null;
	}
	,getKnownFile: function() {
		return this.knownFile;
	}
	,setConstFile: function(content) {
		this.constFile = content;
	}
	,hasConstFile: function() {
		return this.constFile != null;
	}
	,getConstFile: function() {
		return this.constFile;
	}
};
var Std = function() { };
Std.__name__ = true;
Std.string = function(s) {
	return js_Boot.__string_rec(s,"");
};
var StringTools = function() { };
StringTools.__name__ = true;
StringTools.replace = function(s,sub,by) {
	return s.split(sub).join(by);
};
var haxe_Exception = function(message,previous,native) {
	Error.call(this,message);
	this.message = message;
	this.__previousException = previous;
	this.__nativeException = native != null ? native : this;
};
haxe_Exception.__name__ = true;
haxe_Exception.thrown = function(value) {
	if(((value) instanceof haxe_Exception)) {
		return value.get_native();
	} else if(((value) instanceof Error)) {
		return value;
	} else {
		var e = new haxe_ValueException(value);
		return e;
	}
};
haxe_Exception.__super__ = Error;
haxe_Exception.prototype = $extend(Error.prototype,{
	get_native: function() {
		return this.__nativeException;
	}
});
var haxe_ValueException = function(value,previous,native) {
	haxe_Exception.call(this,String(value),previous,native);
	this.value = value;
};
haxe_ValueException.__name__ = true;
haxe_ValueException.__super__ = haxe_Exception;
haxe_ValueException.prototype = $extend(haxe_Exception.prototype,{
});
var haxe_ds_List = function() {
	this.length = 0;
};
haxe_ds_List.__name__ = true;
haxe_ds_List.prototype = {
	add: function(item) {
		var x = new haxe_ds__$List_ListNode(item,null);
		if(this.h == null) {
			this.h = x;
		} else {
			this.q.next = x;
		}
		this.q = x;
		this.length++;
	}
	,join: function(sep) {
		var s_b = "";
		var first = true;
		var l = this.h;
		while(l != null) {
			if(first) {
				first = false;
			} else {
				s_b += sep == null ? "null" : "" + sep;
			}
			s_b += Std.string(l.item);
			l = l.next;
		}
		return s_b;
	}
};
var haxe_ds__$List_ListNode = function(item,next) {
	this.item = item;
	this.next = next;
};
haxe_ds__$List_ListNode.__name__ = true;
var haxe_ds_StringMap = function() {
	this.h = Object.create(null);
};
haxe_ds_StringMap.__name__ = true;
var haxe_iterators_ArrayIterator = function(array) {
	this.current = 0;
	this.array = array;
};
haxe_iterators_ArrayIterator.__name__ = true;
haxe_iterators_ArrayIterator.prototype = {
	hasNext: function() {
		return this.current < this.array.length;
	}
	,next: function() {
		return this.array[this.current++];
	}
};
var js_Boot = function() { };
js_Boot.__name__ = true;
js_Boot.__string_rec = function(o,s) {
	if(o == null) {
		return "null";
	}
	if(s.length >= 5) {
		return "<...>";
	}
	var t = typeof(o);
	if(t == "function" && (o.__name__ || o.__ename__)) {
		t = "object";
	}
	switch(t) {
	case "function":
		return "<function>";
	case "object":
		if(((o) instanceof Array)) {
			var str = "[";
			s += "\t";
			var _g = 0;
			var _g1 = o.length;
			while(_g < _g1) {
				var i = _g++;
				str += (i > 0 ? "," : "") + js_Boot.__string_rec(o[i],s);
			}
			str += "]";
			return str;
		}
		var tostr;
		try {
			tostr = o.toString;
		} catch( _g ) {
			return "???";
		}
		if(tostr != null && tostr != Object.toString && typeof(tostr) == "function") {
			var s2 = o.toString();
			if(s2 != "[object Object]") {
				return s2;
			}
		}
		var str = "{\n";
		s += "\t";
		var hasp = o.hasOwnProperty != null;
		var k = null;
		for( k in o ) {
		if(hasp && !o.hasOwnProperty(k)) {
			continue;
		}
		if(k == "prototype" || k == "__class__" || k == "__super__" || k == "__interfaces__" || k == "__properties__") {
			continue;
		}
		if(str.length != 2) {
			str += ", \n";
		}
		str += s + k + " : " + js_Boot.__string_rec(o[k],s);
		}
		s = s.substring(1);
		str += "\n" + s + "}";
		return str;
	case "string":
		return o;
	default:
		return String(o);
	}
};
if(typeof(performance) != "undefined" ? typeof(performance.now) == "function" : false) {
	HxOverrides.now = performance.now.bind(performance);
}
String.__name__ = true;
Array.__name__ = true;
js_Boot.__toStr = ({ }).toString;
FastaAlignmentParser.authorizedCharacters = (function($this) {
	var $r;
	var _g = new haxe_ds_StringMap();
	_g.h["A"] = true;
	_g.h["T"] = true;
	_g.h["G"] = true;
	_g.h["C"] = true;
	_g.h["N"] = true;
	_g.h["-"] = true;
	_g.h["?"] = true;
	_g.h["R"] = false;
	_g.h["Y"] = false;
	_g.h["M"] = false;
	_g.h["K"] = false;
	_g.h["W"] = false;
	_g.h["S"] = false;
	$r = _g;
	return $r;
}(this));
SeqPhase1.map1 = (function($this) {
	var $r;
	var _g = new haxe_ds_StringMap();
	_g.h["W"] = "A";
	_g.h["S"] = "C";
	_g.h["K"] = "T";
	_g.h["M"] = "A";
	_g.h["Y"] = "C";
	_g.h["R"] = "A";
	$r = _g;
	return $r;
}(this));
SeqPhase1.map2 = (function($this) {
	var $r;
	var _g = new haxe_ds_StringMap();
	_g.h["W"] = "T";
	_g.h["S"] = "G";
	_g.h["K"] = "G";
	_g.h["M"] = "C";
	_g.h["Y"] = "T";
	_g.h["R"] = "G";
	$r = _g;
	return $r;
}(this));
SeqPhase1.code = (function($this) {
	var $r;
	var _g = new haxe_ds_StringMap();
	_g.h["A"] = "1";
	_g.h["C"] = "2";
	_g.h["G"] = "3";
	_g.h["T"] = "4";
	_g.h["?"] = "?";
	_g.h["N"] = "?";
	_g.h["-"] = "0";
	$r = _g;
	return $r;
}(this));
