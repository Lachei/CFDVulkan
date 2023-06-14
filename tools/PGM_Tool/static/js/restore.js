import {set_size} from './modules/functions.mjs';
import {pgm_dict_flipped} from './modules/dicts.mjs';
import { logger } from './modules/logger.mjs';
import { table } from './modules/table.mjs';

if (window.hasOwnProperty( "array_upload" )) {
    if ((rows_upload != 'None') && (cols_upload != 'None')) {
        set_size(rows_upload, cols_upload);
    };

    if (file_type.localeCompare("pgm") == 0) {
        // Iterate through array
        for (let row = 0; row < rows_upload; row++) {
            for (let col = 0; col < cols_upload; col++) {
                let old = array_upload[row][col];
                // Perform replacement of old values with new values
                array_upload[row][col] = pgm_dict_flipped[old].toUpperCase()
            }
            // Set Row Data for each row of upload array
            table.setRowData(row, array_upload[row]);
        }
        logger.log(`Restored PGM File: ${filename}`);
    }
    else if (file_type.localeCompare("csv") == 0) {
        for (let row = 0; row < rows_upload; row++) {
            // Set Row Data for each row of upload array
            table.setRowData(row, array_upload[row]);
        }
        logger.log(`Restored CSV File: ${filename}`);
    }
}
