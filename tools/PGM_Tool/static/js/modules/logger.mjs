/**
 * @author Moritz Gnisia <moritz.gnisia@tum.de>
 *
 * Created at: 2020-01-18
**/

// Logger Idea based on the following script: https://code.sololearn.com/Wn395PKwfxJ4/#html

// Creating Logger Object
export let logger = {};

// Getting div to insert logs
let logger_object = document.getElementById("logger");

// Basic Formating for Logging Output
function format_output(text) {
     /* Time Formating with Padding 
    Resource: https://stackoverflow.com/questions/3313875/javascript-date-ensure-getminutes-gethours-getseconds-puts-0-in-front-i */
    let d = new Date()
    let h = `${d.getHours()}`.padStart(2, '0')
    let m = `${d.getMinutes()}`.padStart(2, '0')
    let s = `${d.getSeconds()}`.padStart(2, '0')

    // Create Output: Combination of time + Text Input
    return h + ":" + m + ":" + s + " - " + text;
};

// This function appends the text output to the div object
function append_text(text, color="") {
    let element = document.createElement("div");
    if (color.length > 2) {
        element.style.color = 'red';
    }
    element.appendChild(text);
    logger_object.appendChild(element);
}

// Adding log method from our logger object
logger.log = text => {
    let txt = document.createTextNode(format_output(text));
    append_text(txt);
}

// Adding log method from our logger object
logger.warn = text => {
    let txt = document.createTextNode(format_output(text));
    append_text(txt,"red");
}