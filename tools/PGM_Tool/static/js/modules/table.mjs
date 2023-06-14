// Table Setup
import { color_dict } from './dicts.mjs';
import { logger } from './logger.mjs';

// Definitions of Functions
var options = {
    minDimensions:[1,1],
    allowExport:true,
    tableOverflow: true,
    tableWidth: "1000px",
    tableHeight: "300px",
    updateTable:function(instance, cell, col, row, val, label, cellName) {
    // Cell Formating
    var text_value = cell.innerText.toLocaleLowerCase();
    if (text_value.length >= 1)  {
        if (text_value in color_dict) {
            cell.innerText=cell.innerText.toLocaleUpperCase();
            cell.style.backgroundColor = color_dict[text_value][0];
            cell.style.color = color_dict[text_value][1];
        }
        else {
            logger.warn(`The value ${text_value} entered for the cell doesn't match with the legend. It will be removed.`);
            cell.innerText='';
            cell.style.backgroundColor = '#ffffff';
        }
    }
    }
}

export var table = jexcel(document.getElementById('spreadsheet'), options);
