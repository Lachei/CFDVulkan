/**
 * @author Moritz Gnisia <moritz.gnisia@tum.de>
 *
 * Created at: 2020-01-18
**/
import { logger } from './modules/logger.mjs';
import { color_dict } from './modules/dicts.mjs';
import { table } from './modules/table.mjs';
import { set_size_user, download_pgm, set_fluid, download_csv} from './modules/functions.mjs';

// For Loop through the color_dict
// This creates the legend on the right hand side for index.html
for (var key in color_dict) {
    if (color_dict.hasOwnProperty(key)) {
        $('#legend').find('tbody:last').append('<tr><td>' + key.toLocaleUpperCase() + "</td><td class='legend' bgcolor='" +
            color_dict[key][0] + "'></td><td>" + color_dict[key][2] + "</td></tr>");
    }
}
logger.log("Legend initialized");

// Initialize Domain
let size = 4;
table.insertRow(size, 0, true);
table.insertColumn(size,0, true);
logger.log("Initialized Domain to a size of 5x5");

// Init Buttons

// Create Button Click function for Set Size
let BTN_set_size = document.querySelector('#set_size');
BTN_set_size.addEventListener('click', () => {
    set_size_user();
});

// Create Button Click function for Fluid Setting
let BTN_set_fluid = document.querySelector('#set_fluid');
BTN_set_fluid.addEventListener('click', () => {
    set_fluid();
});

// Create Button Click function for Download PGM
let BTN_download_pgm = document.querySelector('#download_pgm');
BTN_download_pgm.addEventListener('click', () => {
    download_pgm();
});

// Create Button Click function for Download CSV File
let BTN_download_csv = document.querySelector('#download_csv');
BTN_download_csv.addEventListener('click', () => {
    download_csv();
});