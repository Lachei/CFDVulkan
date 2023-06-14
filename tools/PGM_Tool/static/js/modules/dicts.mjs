// Creating a dictionary for mapping of the cell types to the colors
export let color_dict = {
    "f": ["#0097e6", 'black', 'Fluid cell'],
    "o": ["#718093", 'white', 'Obstacle cell'],
    "i": ["#dcdde1", 'black', 'Inlet cell'],
    "ot": ["#130f40", 'white', 'Outlet cell'],
};

export let pgm_dict = {
    "f":  1,
    "o":  2,
    "i":  3,
    "ot": 4,
};

export let pgm_dict_flipped = swap(pgm_dict);

// Swap Function from the following SOF entry: https://stackoverflow.com/questions/23013573/swap-key-with-value-json/23013744
function swap(json){
    var ret = {};
    for(var key in json){
      ret[json[key]] = key;
    }
    return ret;
}
  