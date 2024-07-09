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
class Step2
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

    public static function setConstFileContent(str) {
        var displayStr:String = "" + str;
        if(str.length > 15) {
            displayStr = str.substr(0, 12) + "...";
        }
        Browser.document.getElementById("passedSpan").style.display = "block";
        var cell = Browser.document.getElementById("passedConstFileContent");
        cell.textContent = displayStr;
        var changeLink = Browser.document.getElementById("changePassedConstFile");
        changeLink.onclick = resetConstFileUpload;
        Browser.document.getElementById("passedDescription").style.display = "block";
        var constFile = Browser.document.getElementById("constFile");
        constFile.style.display = "none";
        var textField:js.html.InputElement = cast Browser.document.getElementById("completePassedConstFileContent");
        textField.value = str;
    }

    public static function clearURI() {
        Browser.window.history.replaceState({}, Browser.document.title, Browser.window.location.pathname);
    }

    public static function resetForm() {
        var textField:js.html.InputElement = cast Browser.document.getElementById("phaseOutFile");
        textField.value = "";
        resetConstFileUpload();
        var textField2:js.html.InputElement = cast Browser.document.getElementById("constFile");
        textField2.value = "";
        var checkbox:js.html.InputElement = cast Browser.document.getElementById("sortSequences");
        checkbox.checked = true;
        var checkbox2:js.html.InputElement = cast Browser.document.getElementById("reduceSequences");
        checkbox2.checked = true;
    }

    public static function handleResults(result) {
        Browser.document.getElementById("resultArea").style.display = "block";
        var b64 = Browser.window.btoa(result);
        var link:js.html.LinkElement = cast Browser.document.getElementById("downloadLink");
        link.href = 'data:text/plain;base64,\n'+b64;
    }

    public static function run(sort, reduce, outFile, constFile) {
        var result = SeqPhase2.parse(outFile, constFile).getFasta(reduce, sort);
        Browser.window.setTimeout(function() {
            handleResults(result);
        }, 0);
    }

    public static function runReadConstFile(sort, reduce, outFileContent) {
        var constFile:js.html.InputElement = cast Browser.document.getElementById("constFile");
        var passedConstFileContentField:js.html.InputElement = cast Browser.document.getElementById("completePassedConstFileContent");
        if(passedConstFileContentField.value != "") {
            var constFileContentField:js.html.InputElement = cast Browser.document.getElementById("completePassedConstFileContent");
            var constFileContent = constFileContentField.value;
            Browser.window.setTimeout(function() {
                run(sort, reduce, outFileContent, constFileContent);
            }, 0);
        } else if(constFile.files.length == 0) {
            Browser.window.setTimeout(function() {
                run(sort, reduce, outFileContent, "");
            }, 0);
        } else if(constFile.files.length == 1) {
            var reader = new FileReader();
            reader.onload = function(data) {
                var constFileContent = data.target.result;
                Browser.window.setTimeout(function() {
                    run(sort, reduce, outFileContent, constFileContent);
                }, 0);
            }
            reader.readAsText(constFile.files[0], "ISO-8859-1");
        }
    }

    public static function runReadOutFile(sort, reduce) {
        var outFileInput:js.html.InputElement = cast Browser.document.getElementById("phaseOutFile");
        if(outFileInput.files.length == 1) {
            var reader = new FileReader();
            reader.onload = function(data) {
                var outFileContent = data.target.result;
                Browser.window.setTimeout(function() {
                    runReadConstFile(sort, reduce, outFileContent);
                }, 0);
            }
            reader.readAsText(outFileInput.files[0], "ISO-8859-1");
        } else {
            Browser.window.alert("Missing .out / .out_pairs file");
        }
    }

    public static function runJob() {
        var sortE:js.html.InputElement = cast Browser.document.getElementById("sortSequences");
        var sort = sortE.checked ;
        var reduceE:js.html.InputElement = cast Browser.document.getElementById("reduceSequences");
        var reduce = reduceE.checked;
        runReadOutFile(sort, reduce);
    }

    public static function main() {
        var uri = Browser.window.location.search;
        if(uri != null && uri != "") {
            if(StringTools.startsWith(uri, "?")) {
                uri = uri.substr(1);
            }
            var parts = uri.split("&");
            for(i in 0...parts.length) {
                if(StringTools.startsWith(parts[i], "constFileContent=")) {
                    var constFileContent = StringTools.trim(parts[i].substr(17)).toUpperCase();
                    if(constFileContent != null && constFileContent != "") {
                        setConstFileContent(constFileContent);
                        clearURI();
                    }
                }
            }
        }
        Browser.document.getElementById("resetJob2Button").onclick = resetForm;
        Browser.document.getElementById("runJob2Button").onclick = runJob;
    }
}
