// Function Definitions
import { logger } from './logger.mjs';
import { table } from './table.mjs';
import { pgm_dict } from './dicts.mjs';

export function set_size(rows, cols) {
    let size = get_size();
    let rows_current = size[0];
    let cols_current = size[1];
   
    if (rows_current > rows) {
        table.deleteRow(0, rows_current - rows);
    }
    else if (rows_current < rows) {
        table.insertRow(rows - rows_current, 0, true);
    }
    if (cols_current > cols) {
        table.deleteColumn(0, cols_current - cols);
    }
    else if (cols_current < cols) {
        table.insertColumn(cols - cols_current,0, true);
    }
    logger.log(`Set new size: Rows ${rows} and  Columns ${cols}`);
}

export function set_size_user() {
    let rows = $('#rows').val();
    let cols = $('#cols').val();
    set_size(rows, cols);
}

export function download_pgm() {
    /* Get Table data for all cells, therefore boolean set to false
    Otherwise only highlighted cells would be return
    */
    let size = get_size();
    let rows = size[0];
    let cols = size[1];
    
    let header = 'P2 \n' + `${cols} ${rows} \n`;
    let array_str = "";
    for (let row = 0; row < rows; row++) {
        for (let col = 0; col < cols; col++) {
            let val = table.getValueFromCoords(row, col).toLowerCase();
            if (val in pgm_dict) {
                array_str += pgm_dict[val] + " ";
            }
            else {
                array_str += " " + val + " ";
            }
        }
        array_str += '\n';
    }
    // Merge Header and array
    let pgm = header + array_str;

    // Date Configuration for Filename
    let date = new Date();
    let d = `${date.getDate()}`.padStart(2, '0');
    let m = `${date.getMonth()+1}`.padStart(2, '0'); // Date starts with 0 till 11, therefore +1
    let y = `${date.getFullYear()}`;

    // Compose Filename, e.g. 2020-01-19_Config.pgm
    let filename = y + "-" + m + "-" + d + "_" + "Config.pgm";
    
    // Prepare and start download
    let download = new Blob([pgm], {type: "text/plain;charset=utf-8"});
    logger.log("Downloading Config as PGM File");
    saveAs(download, filename);
}

export function download_csv(){
    table.download();
    logger.log("Downloading Config as CSV File");
}

function writedata() {
    let data = table.getData(false);
    console.log(data);
}

function get_size() {
    let data = table.getData();
    let rows = data.length;
    let cols = data[0].length;
    return [rows, cols];
}

export function set_fluid() {
    var size = get_size();
    var cols = size[1];
    for(var i=0; i < cols; i++) {
        table.setColumnData(i, table.getColumnData(i).fill('F'));
    }
    logger.log("Initialized all cells to Fluid cells");
}

