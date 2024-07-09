/**
 * Copyright (c) 2024, Yann Spöri
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import haxe.ds.Vector;
import haxe.Exception;

import js.Browser;
import js.html.FileReader;
import js.html.DOMElement;

/**
 * Some Haxe code (will be compiled to JavaScript) for step1.
 *
 * @author Yann Spoeri
 */
class Step1
{
    public static function errFunction(msg, url, line:Int, col, error):Dynamic {
        Browser.window.alert("Error " + msg + ", line " + line);
        return null;
    }
    
    public static function resetConstFileUpload() {
        Browser.document.getElementById("passedSpan").style.display = "none";
        var cell = Browser.document.getElementById("passedConstFileContent");
        cell.textContent = "";
        var changeLink = Browser.document.getElementById("changePassedConstFile");
        changeLink.onclick = resetConstFileUpload;
        Browser.document.getElementById("passedDescription").style.display = "none";
        var constFile = Browser.document.getElementById("constFile");
        constFile.style.display = "block";
        var textField:js.html.InputElement = cast Browser.document.getElementById("completePassedConstFileContent");
        textField.value = "";
    }
    
    public static function resetForm() {
        var textField1:js.html.InputElement = cast Browser.document.getElementById("alig1");
        textField1.value = "";
        var textField2:js.html.InputElement = cast Browser.document.getElementById("alig2");
        textField2.value = "";
        var textField3:js.html.InputElement = cast Browser.document.getElementById("alig3");
        textField3.value = "";
    }
    
    public static function copySuggestedCommand() {
        var copyText:js.html.InputElement = cast Browser.document.getElementById('suggestedCommandField');
        copyText.select();
        Browser.document.execCommand("copy");
    }
    
    public static function resetResultArea() {
        for(i in 1...4) {
            Browser.document.getElementById("align" + i + "FileErrorsArea").style.display = "none";
            Browser.document.getElementById("align" + i + "FileErrorsList").innerHTML = "";
            Browser.document.getElementById("align" + i + "FileWarningsArea").style.display = "none";
            Browser.document.getElementById("align" + i + "FileWarningsList").innerHTML = "";
        }
        Browser.document.getElementById("generalErrorsArea").style.display = "none";
        Browser.document.getElementById("generalErrorsList").innerHTML = "";
        Browser.document.getElementById("generalWarningsArea").style.display = "none";
        Browser.document.getElementById("generalWarningsList").innerHTML = "";
        Browser.document.getElementById("notesArea").innerHTML = "";
        Browser.document.getElementById("outFileResult").style.display = "none";
        Browser.document.getElementById("knownFileResult").style.display = "none";
        Browser.document.getElementById("constFileResult").style.display = "none";
        Browser.document.getElementById("suggestedCommandArea").style.display = "none";
        Browser.document.getElementById("proceedArea").style.display = "none";
    }

    public static function handleResults(result:SeqPhase1Result) {
        resetResultArea();
        for(k in 1...4) {
            if(result.hasAlignErrors(k)) {
                Browser.document.getElementById("align" + k + "FileErrorsArea").style.display = "block";
                var errors = result.getAlignErrors(k);
                for(i in 0...errors.length) {
                    var ele = Browser.document.createElement("li");
                    ele.innerText = errors[i];
                    Browser.document.getElementById("align" + k + "FileErrorsList").append(ele);
                }
            }
            if(result.hasAlignWarn(k)) {
                Browser.document.getElementById("align" + k + "FileWarningsArea").style.display = "block";
                var warnings = result.getAlignWarn(k);
                for(i in 0...warnings.length) {
                    var ele = Browser.document.createElement("li");
                    ele.innerText = warnings[i];
                    Browser.document.getElementById("align" + k + "FileWarningsList").append(ele);
                }
            }
        }
        if(result.hasGeneralErrors()) {
            Browser.document.getElementById("generalErrorsArea").style.display = "block";
            var errors = result.getGeneralErrors();
            for(i in 0...errors.length) {
                var ele = Browser.document.createElement("li");
                ele.innerText = errors[i];
                Browser.document.getElementById("generalErrorsList").append(ele);
            }
        }
        if(result.hasGeneralWarn()) {
            Browser.document.getElementById("generalWarningsArea").style.display = "block";
            var warnings = result.getGeneralWarn();
            for(i in 0...warnings.length) {
                var ele = Browser.document.createElement("li");
                ele.innerText = warnings[i];
                Browser.document.getElementById("generalWarningsList").append(ele);
            }
        }
        if(result.hasNotes()) {
            Browser.document.getElementById("notesArea").style.display = "block";
            var notes = result.getNotes();
            for(i in 0...notes.length) {
                var ele = Browser.document.createElement("br");
                Browser.document.getElementById("notesArea").append(ele);
                var ele = Browser.document.createElement("span");
                ele.innerText = notes[i];
                Browser.document.getElementById("notesArea").append(ele);
            }
        }
        if(result.hasSuggestedCommand()) {
            Browser.document.getElementById("suggestedCommandArea").style.display = "block";
            var textField1:js.html.InputElement = cast Browser.document.getElementById("suggestedCommandField");
            textField1.value = result.getSuggestedPhaseCommand();
        }
        if(result.hasInpFile()) {
            Browser.document.getElementById("outFileResult").style.display = "block";
            var link:js.html.LinkElement = cast Browser.document.getElementById("downloadInpFileLink");
            var b64 = Browser.window.btoa(result.getInpFile());
            link.href = 'data:text/plain;base64,\n'+b64;
        }
        if(result.hasKnownFile()) {
            Browser.document.getElementById("knownFileResult").style.display = "block";
            var link:js.html.LinkElement = cast Browser.document.getElementById("downloadKnownFileLink");
            var b64 = Browser.window.btoa(result.getKnownFile());
            link.href = 'data:text/plain;base64,\n'+b64;
        }
        if(result.hasConstFile()) {
            Browser.document.getElementById("constFileResult").style.display = "block";
            var link:js.html.LinkElement = cast Browser.document.getElementById("downloadConstFileLink");
            var b64 = Browser.window.btoa(result.getConstFile());
            link.href = 'data:text/plain;base64,\n'+b64;
        }
        if(result.hasInpFile()) {
            Browser.document.getElementById("proceedArea").style.display = "block";
            if(result.hasConstFile()) {
                var ele:js.html.LinkElement = cast Browser.document.getElementById("finalLink");
                ele.href = "step2.html?constFileContent=" + result.getConstFile();
            } else {
                var ele:js.html.LinkElement = cast Browser.document.getElementById("finalLink");
                ele.href = "step2.html";
            }
        }
    }
    
    public static function run(align1Content, align2Content, align3Content) {
        var result = SeqPhase1.doIt(align1Content, align2Content, align3Content);
        Browser.window.setTimeout(function() {
            handleResults(result);
        }, 0);
    }
    
    public static function runReadAlign3(align1Content, align2Content) {
        var fObj:js.html.InputElement = cast Browser.document.getElementById("alig3");
        if(fObj.files.length == 1) {
            var reader = new FileReader();
            reader.onload = function(data) {
                var fileContent = data.target.result;
                Browser.window.setTimeout(function() {
                    run(align1Content, align2Content, fileContent);
                }, 0);
            }
            reader.readAsText(fObj.files[0], "ISO-8859-1");
        } else {
            Browser.window.setTimeout(function() {
                run(align1Content, align2Content, null);
            }, 0);
        }    
    }
    
    public static function runReadAlign2(align1Content) {
        var fObj:js.html.InputElement = cast Browser.document.getElementById("alig2");
        if(fObj.files.length == 1) {
            var reader = new FileReader();
            reader.onload = function(data) {
                var fileContent = data.target.result;
                Browser.window.setTimeout(function() {
                    runReadAlign3(align1Content, fileContent);
                }, 0);
            }
            reader.readAsText(fObj.files[0], "ISO-8859-1");
        } else {
            if(align1Content == null) {
                Browser.window.alert("Missing input data! Either an alignment of sequences from homozygous individuals and from heterozygotes to be phased or an alignment of fake haplotype pairs from heterozygotes to be phased has to be provided. You can also provide both files!");
            } else {
                Browser.window.setTimeout(function() {
                    runReadAlign3(align1Content, null);
                }, 0);
            }      
        }    
    }

    public static function runReadAlign1() {
        var fObj:js.html.InputElement = cast Browser.document.getElementById("alig1");
        if(fObj.files.length == 1) {
            var reader = new FileReader();
            reader.onload = function(data) {
                var fileContent = data.target.result;
                Browser.window.setTimeout(function() {
                    runReadAlign2(fileContent);
                }, 0);
            }
            reader.readAsText(fObj.files[0], "ISO-8859-1");
        } else {
            Browser.window.setTimeout(function() {
                runReadAlign2(null);
            }, 0);        
        }    
    }

    public static function main() {
        Browser.window.onerror = errFunction;
        Browser.document.getElementById("resetJob1Button").onclick = resetForm;
        Browser.document.getElementById("copySuggestedCommand").onclick = copySuggestedCommand;
        Browser.document.getElementById("submitJob1Button").onclick = runReadAlign1;
    }
}
