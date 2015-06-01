<?php
class M
{
    public $time;
    public $pressure;
    public $t1;
    public $t2;
    public $t3;
    public $t4;
    public $t5;
    public $t6;
    public $t_ext;
    public $m1;
    public $m2;
    public $m3;
}

$measurements = array();

function get_measurements($file)
{
    $content = file_get_contents($file);
    
    $output = array();
    preg_match_all("#([:digit]*)/([:digit]*)/([:digit]*)(?:.*)([:digit]*):([:digit]*):([:digit]*)(?:.*)(-?)([:digit]*)(?:.*)(-?)([:digit]*)(?:.*)(-?)([:digit]*)(?:.*)(-?)([:digit]*)(?:.*)(-?)([:digit]*)(?:.*)(-?)([:digit]*)(?:.*)(-?)([:digit]*)(?:.*)(-?)([:digit]*)#sU", $content, $output);
    
    return $measurements;
}

echo get_measurements("IBL_2015_05_29_14_12_02.txt");
?>
