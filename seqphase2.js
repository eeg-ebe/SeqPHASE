function $extend(from, fields) {
	var proto = Object.create(from);
	for (var name in fields) proto[name] = fields[name];
	if( fields.toString !== Object.prototype.toString ) proto.toString = fields.toString;
	return proto;
}
var EReg = function(r,opt) {
	this.r = new RegExp(r,opt.split("u").join(""));
};
EReg.__name__ = true;
EReg.prototype = {
	split: function(s) {
		var d = "#__delim__#";
		return s.replace(this.r,d).split(d);
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
var Individual = function(name,seq1,seq2,index,p) {
	this.name = name;
	this.seq1 = seq1;
	this.seq2 = seq2;
	this.index = index;
	this.p = p;
};
Individual.__name__ = true;
Individual.prototype = {
	getName: function() {
		return this.name;
	}
	,getSequence1: function() {
		return this.seq1;
	}
	,getSequence2: function() {
		return this.seq2;
	}
	,getIndex: function() {
		return this.index;
	}
	,getProbability: function() {
		return this.p;
	}
	,isHomozygous: function() {
		return this.seq1 == this.seq2;
	}
	,isHeterozygous: function() {
		return this.seq1 != this.seq2;
	}
	,getFasta: function(reduce,outFile) {
		var result = new haxe_ds_List();
		if(this.p == 1.0) {
			if(this.isHomozygous() && reduce) {
				result.add(">" + this.name);
				result.add(this.seq1);
			} else {
				result.add(">" + this.name + "a");
				result.add(this.seq1);
				result.add(">" + this.name + "b");
				result.add(this.seq2);
			}
		} else if(this.isHomozygous() && reduce) {
			result.add(">" + this.name + "-" + this.index + "(" + this.p + ")");
			result.add(this.seq1);
		} else {
			result.add(">" + this.name + "-" + this.index + "a(" + this.p + ")");
			result.add(this.seq1);
			result.add(">" + this.name + "-" + this.index + "b(" + this.p + ")");
			result.add(this.seq2);
		}
		return result.join("\n");
	}
};
var Individuals = function() {
	this.inds = [];
};
Individuals.__name__ = true;
Individuals.prototype = {
	addIndividual: function(ind) {
		this.inds.push(ind);
	}
	,setOutFile: function(b) {
		this.outFile = b;
	}
	,getFasta: function(reduce,sort) {
		var result = new haxe_ds_List();
		this.inds.sort(function(i1,i2) {
			if(i1.getName() < i2.getName()) {
				return -1;
			}
			if(i1.getName() > i2.getName()) {
				return 1;
			}
			return i1.getIndex() - i2.getIndex();
		});
		var _g = 0;
		var _g1 = this.inds;
		while(_g < _g1.length) {
			var ind = _g1[_g];
			++_g;
			result.add(ind.getFasta(reduce,this.outFile));
		}
		return result.join("\n");
	}
};
Math.__name__ = true;
var SeqPhase2 = function() { };
SeqPhase2.__name__ = true;
SeqPhase2.removePossibleLineEndingR = function(line) {
	if(line.length >= 1) {
		if(line.charAt(0) == "\r") {
			line = HxOverrides.substr(line,1,null);
		}
		if(line.charAt(line.length - 1) == "\r") {
			line = HxOverrides.substr(line,0,line.length - 1);
		}
	}
	return StringTools.trim(line);
};
SeqPhase2.genConstString = function(length) {
	var result = new haxe_ds_List();
	var _g = 0;
	var _g1 = length;
	while(_g < _g1) {
		var i = _g++;
		result.add(".");
	}
	return result.join("");
};
SeqPhase2.nSplit = function(s) {
	var result = new haxe_ds_List();
	var _g = 0;
	var _g1 = s.length;
	while(_g < _g1) {
		var i = _g++;
		var char = s.charAt(i);
		if(char == " ") {
			continue;
		}
		result.add(char);
	}
	return result;
};
SeqPhase2.nrsToSeq = function(nrs,$const) {
	var nrs1 = SeqPhase2.nSplit(nrs);
	if($const == null || $const == "") {
		$const = SeqPhase2.genConstString(nrs1.length);
	}
	var codeS = new haxe_ds_List();
	var _g_head = nrs1.h;
	while(_g_head != null) {
		var val = _g_head.item;
		_g_head = _g_head.next;
		var part = val;
		var i = Std.parseInt(part);
		var r = SeqPhase2.code.h[i];
		if(r == null) {
			throw haxe_Exception.thrown("Error: The .out/.out_pairs file contains an allele not present in the input dataset (" + i + ", " + Std.string(nrs1) + ")");
		}
		codeS.add(r);
	}
	var doesNotMatch = false;
	var result = new haxe_ds_List();
	var _g = 0;
	var _g1 = $const.length;
	while(_g < _g1) {
		var i = _g++;
		var chr = $const.charAt(i);
		if(chr == ".") {
			var nextCodeS = codeS.pop();
			if(nextCodeS == null) {
				doesNotMatch = true;
				break;
			}
			result.add(nextCodeS);
		} else {
			result.add(chr);
		}
	}
	if(doesNotMatch || codeS.length != 0) {
		throw haxe_Exception.thrown("Error: The data in the const and in the input file do not match; please check input data. (" + (doesNotMatch == null ? "null" : "" + doesNotMatch) + ", " + codeS.length + ")");
	}
	return result.join("");
};
SeqPhase2.parsePairFile = function(fileContent,$const) {
	if($const == null) {
		$const = "";
	}
	var result = new Individuals();
	result.setOutFile(false);
	var lineNo = 0;
	var currentIndividualName = "";
	var index = 0;
	var _g = 0;
	var _g1 = fileContent.split("\n");
	while(_g < _g1.length) {
		var line = _g1[_g];
		++_g;
		line = SeqPhase2.removePossibleLineEndingR(line);
		++lineNo;
		if(line != "") {
			if(HxOverrides.substr(line,0,5) == "IND: ") {
				currentIndividualName = HxOverrides.substr(line,5,null);
				index = 1;
			} else {
				var parts = line.split(" , ");
				if(parts.length != 3) {
					throw haxe_Exception.thrown("Please check line " + lineNo + " in the input file!");
				}
				var seq1 = SeqPhase2.nrsToSeq(parts[0],$const);
				var seq2 = SeqPhase2.nrsToSeq(parts[1],$const);
				var probability = parseFloat(parts[2]);
				var ind = new Individual(currentIndividualName,seq1,seq2,index++,probability);
				result.addIndividual(ind);
			}
		}
	}
	return result;
};
SeqPhase2.stripBraces = function(c) {
	if(c.charAt(0) == "(") {
		c = HxOverrides.substr(c,1,null);
	}
	if(c.charAt(c.length - 1) == ")") {
		c = HxOverrides.substr(c,0,c.length - 1);
	}
	return c;
};
SeqPhase2.getSumCode = function(c1,c2) {
	c1 = SeqPhase2.stripBraces(c1);
	c2 = SeqPhase2.stripBraces(c2);
	var code1 = Std.parseInt(c1);
	var code2 = Std.parseInt(c2);
	var code1S = SeqPhase2.code.h[code1];
	var code2S = SeqPhase2.code.h[code2];
	var codeSum = SeqPhase2.supercode.h[code1S] + SeqPhase2.supercode.h[code2S];
	return SeqPhase2.revSupercode.h[codeSum];
};
SeqPhase2.processLines = function(indName,seqLine1,seqLine2,$const,lineNo) {
	var result = new haxe_ds_List();
	var r = new EReg(" +","g");
	var pLine1 = r.split(seqLine1);
	var pLine2 = r.split(seqLine2);
	if(pLine1.length != pLine2.length) {
		throw haxe_Exception.thrown("Error: Sequences have different lengths. Please check sequences of individual " + indName + " (line " + (lineNo - 2) + "ff.) in input file.");
	}
	if($const == null || $const == "") {
		$const = SeqPhase2.genConstString(pLine1.length);
	}
	var r1 = new haxe_ds_List();
	var r2 = new haxe_ds_List();
	var _g = 0;
	var _g1 = pLine1.length;
	while(_g < _g1) {
		var i = _g++;
		var c1 = pLine1[i];
		var c2 = pLine2[i];
		var processed1 = false;
		var processed2 = false;
		if(c1.charAt(0) == "[") {
			r1.add("N");
			processed1 = true;
		}
		if(c2.charAt(0) == "[") {
			r2.add("N");
			processed2 = true;
		}
		if(c1.charAt(0) == "(") {
			var toAdd = SeqPhase2.getSumCode(c1,c2);
			if(toAdd == null || toAdd == "") {
				throw haxe_Exception.thrown("Unexpected input: Parsing error for " + indName + " " + toAdd + " " + i + " " + c1 + " " + c2 + " (1)");
			}
			r1.add(toAdd);
			processed1 = true;
		}
		if(c2.charAt(0) == "(") {
			var toAdd1 = SeqPhase2.getSumCode(c2,c1);
			if(toAdd1 == null || toAdd1 == "") {
				throw haxe_Exception.thrown("Unexpected input: Parsing error for " + indName + " " + toAdd1 + " " + i + " " + c1 + " " + c2 + " (2)");
			}
			r2.add(toAdd1);
			processed2 = true;
		}
		if(!processed1) {
			var toAdd2 = SeqPhase2.code.h[Std.parseInt(c1)];
			if(toAdd2 == null || toAdd2 == "") {
				throw haxe_Exception.thrown("Unexpected input: Parsing error for " + indName + " " + toAdd2 + " " + i + " " + c1 + " " + c2 + " (3)");
			}
			r1.add(toAdd2);
		}
		if(!processed2) {
			var toAdd3 = SeqPhase2.code.h[Std.parseInt(c2)];
			if(toAdd3 == null || toAdd3 == "") {
				throw haxe_Exception.thrown("Unexpected input: Parsing error for " + indName + " " + toAdd3 + " " + i + " " + c1 + " " + c2 + " (4)");
			}
			r2.add(toAdd3);
		}
	}
	var seq1 = new haxe_ds_List();
	var seq2 = new haxe_ds_List();
	var doesNotMatch = false;
	var _g = 0;
	var _g1 = $const.length;
	while(_g < _g1) {
		var i = _g++;
		var chr = $const.charAt(i);
		if(chr == ".") {
			var r1S = r1.pop();
			var r2S = r2.pop();
			if(r1S == null || r2S == null) {
				doesNotMatch = true;
				break;
			}
			seq1.add(r1S);
			seq2.add(r2S);
		} else {
			seq1.add(chr);
			seq2.add(chr);
		}
	}
	if(doesNotMatch || r1.length != 0 || r2.length != 0) {
		throw haxe_Exception.thrown("Error: The data in the const and in the input file do not match; please check input data. (" + (doesNotMatch == null ? "null" : "" + doesNotMatch) + ", " + r1.length + ", " + r2.length + ")");
	}
	var ind = new Individual(indName,seq1.join(""),seq2.join(""),-1,1.0);
	result.add(ind);
	return result;
};
SeqPhase2.parseOutFile = function(fileContent,$const) {
	if($const == null) {
		$const = "";
	}
	var result = new Individuals();
	result.setOutFile(true);
	var lineNo = 0;
	var currentIndividualName = null;
	var currentSeq1 = null;
	var currentSeq2 = null;
	var index = 0;
	var seenBeginBestPairs1 = false;
	var _g = 0;
	var _g1 = fileContent.split("\n");
	while(_g < _g1.length) {
		var line = _g1[_g];
		++_g;
		line = SeqPhase2.removePossibleLineEndingR(line);
		++lineNo;
		if(line == "BEGIN BESTPAIRS1" || line == "BEGIN GENOTYPES") {
			seenBeginBestPairs1 = true;
		} else if(line == "END BESTPAIRS1" || line == "END GENOTYPES") {
			break;
		} else if(seenBeginBestPairs1) {
			var isThatAnEmptyLine = StringTools.trim(line);
			if(isThatAnEmptyLine == "") {
				continue;
			}
			if(currentIndividualName == null) {
				currentIndividualName = line.split(" ")[1];
				if(currentIndividualName == null) {
					currentIndividualName = line;
				}
			} else if(currentSeq1 == null) {
				currentSeq1 = line;
			} else {
				currentSeq2 = line;
				var _g_head = SeqPhase2.processLines(currentIndividualName,currentSeq1,currentSeq2,$const,lineNo).h;
				while(_g_head != null) {
					var val = _g_head.item;
					_g_head = _g_head.next;
					var ind = val;
					result.addIndividual(ind);
				}
				currentIndividualName = null;
				currentSeq1 = null;
			}
		}
	}
	return result;
};
SeqPhase2.parse = function(fileContent,$const) {
	if($const == null) {
		$const = "";
	}
	if(HxOverrides.substr(fileContent,0,10) == "**********") {
		return SeqPhase2.parseOutFile(fileContent,$const);
	} else if(HxOverrides.substr(fileContent,0,5) == "IND: ") {
		return SeqPhase2.parsePairFile(fileContent,$const);
	}
	throw haxe_Exception.thrown("Input file not recognized (this is neither a .out nor a .out_pairs file generated by PHASE)");
};
SeqPhase2.main = function() {
};
var Std = function() { };
Std.__name__ = true;
Std.string = function(s) {
	return js_Boot.__string_rec(s,"");
};
Std.parseInt = function(x) {
	if(x != null) {
		var _g = 0;
		var _g1 = x.length;
		while(_g < _g1) {
			var i = _g++;
			var c = x.charCodeAt(i);
			if(c <= 8 || c >= 14 && c != 32 && c != 45) {
				var nc = x.charCodeAt(i + 1);
				var v = parseInt(x,nc == 120 || nc == 88 ? 16 : 10);
				if(isNaN(v)) {
					return null;
				} else {
					return v;
				}
			}
		}
	}
	return null;
};
var StringTools = function() { };
StringTools.__name__ = true;
StringTools.isSpace = function(s,pos) {
	var c = HxOverrides.cca(s,pos);
	if(!(c > 8 && c < 14)) {
		return c == 32;
	} else {
		return true;
	}
};
StringTools.ltrim = function(s) {
	var l = s.length;
	var r = 0;
	while(r < l && StringTools.isSpace(s,r)) ++r;
	if(r > 0) {
		return HxOverrides.substr(s,r,l - r);
	} else {
		return s;
	}
};
StringTools.rtrim = function(s) {
	var l = s.length;
	var r = 0;
	while(r < l && StringTools.isSpace(s,l - r - 1)) ++r;
	if(r > 0) {
		return HxOverrides.substr(s,0,l - r);
	} else {
		return s;
	}
};
StringTools.trim = function(s) {
	return StringTools.ltrim(StringTools.rtrim(s));
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
var haxe_ds_IntMap = function() {
	this.h = { };
};
haxe_ds_IntMap.__name__ = true;
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
	,pop: function() {
		if(this.h == null) {
			return null;
		}
		var x = this.h.item;
		this.h = this.h.next;
		if(this.h == null) {
			this.q = null;
		}
		this.length--;
		return x;
	}
	,toString: function() {
		var s_b = "";
		var first = true;
		var l = this.h;
		s_b += "{";
		while(l != null) {
			if(first) {
				first = false;
			} else {
				s_b += ", ";
			}
			s_b += Std.string(Std.string(l.item));
			l = l.next;
		}
		s_b += "}";
		return s_b;
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
SeqPhase2.code = (function($this) {
	var $r;
	var _g = new haxe_ds_IntMap();
	_g.h[1] = "A";
	_g.h[2] = "C";
	_g.h[3] = "G";
	_g.h[4] = "T";
	_g.h[-1] = "N";
	_g.h[0] = "-";
	$r = _g;
	return $r;
}(this));
SeqPhase2.supercode = (function($this) {
	var $r;
	var _g = new haxe_ds_StringMap();
	_g.h["A"] = 1;
	_g.h["T"] = 2;
	_g.h["G"] = 4;
	_g.h["C"] = 8;
	_g.h["R"] = 5;
	_g.h["Y"] = 10;
	_g.h["M"] = 9;
	_g.h["K"] = 6;
	_g.h["W"] = 3;
	_g.h["S"] = 12;
	_g.h["V"] = 13;
	_g.h["B"] = 14;
	_g.h["H"] = 11;
	_g.h["D"] = 7;
	_g.h["N"] = 15;
	$r = _g;
	return $r;
}(this));
SeqPhase2.revSupercode = (function($this) {
	var $r;
	var _g = new haxe_ds_IntMap();
	_g.h[1] = "A";
	_g.h[2] = "T";
	_g.h[4] = "G";
	_g.h[8] = "C";
	_g.h[5] = "R";
	_g.h[10] = "Y";
	_g.h[9] = "M";
	_g.h[3] = "W";
	_g.h[12] = "S";
	_g.h[6] = "K";
	_g.h[13] = "V";
	_g.h[14] = "B";
	_g.h[11] = "H";
	_g.h[7] = "D";
	_g.h[15] = "N";
	$r = _g;
	return $r;
}(this));
