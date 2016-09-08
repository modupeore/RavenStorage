<?php
//Very important!!!
require_once('atlas_fns.php');
$base_path = "/home/modupe/public_html/atlas";
$date = shell_exec("date +%Y-%m-%d-%T");
$explodedate = substr($date,0,-1);
?>
<?php //Connect to Database 
function db_connect($host, $db='information_schema') {
  if ($db) {
    $db_conn = new mysqli($host, 'frnakenstein', 'maryshelley', $db);
  } else {
    $db_conn = new mysqli($host, 'frnakenstein', 'maryshelley');
  }

  if (mysqli_connect_errno()) {
    throw new Exception("Connection to database failed: " . mysqli_connect_error());
    return false;
  }
  return $db_conn;
}
?>

<?php //Metadata Page function
function meta_display($action, $result, $primary_key) {
  $num_rows = $result->num_rows;
  echo '<br><table class="metadata"><tr>';
  echo '<th align="left" width=40pt bgcolor="white"><font size="2" color="red">Select All</font><input type="checkbox" id="selectall" onClick="selectAll(this)" /></th>';
  $meta = $result->fetch_field_direct(0); echo '<th class="metadata" id="' . $meta->name . '">' . library_id . '</th>';
  $meta = $result->fetch_field_direct(1); echo '<th class="metadata" id="' . $meta->name . '">' . bird_id . '</th>';
  $meta = $result->fetch_field_direct(2); echo '<th class="metadata" id="' . $meta->name . '">' . species . '</th>';
  $meta = $result->fetch_field_direct(3); echo '<th class="metadata" id="' . $meta->name . '">' . line . '</th>';
  $meta = $result->fetch_field_direct(4); echo '<th class="metadata" id="' . $meta->name . '">' . tissue . '</th>';
  $meta = $result->fetch_field_direct(5); echo '<th class="metadata" id="' . $meta->name . '">' . method . '</th>';
  $meta = $result->fetch_field_direct(6); echo '<th class="metadata" id="' . $meta->name . '">' . index . '</th>';
  $meta = $result->fetch_field_direct(7); echo '<th class="metadata" id="' . $meta->name . '">' . chip_result . '</th>';
  $meta = $result->fetch_field_direct(8); echo '<th class="metadata" id="' . $meta->name . '">' . scientist . '</th>';
  $meta = $result->fetch_field_direct(9); echo '<th class="metadata" id="' . $meta->name . '">' . date . '</th>';
  $meta = $result->fetch_field_direct(10); echo '<th class="metadata" id="' . $meta->name . '">' . notes . '</th>';
  $meta = $result->fetch_field_direct(11); echo '<th class="metadata" id="' . $meta->name . '">' . status . '</th></tr>';

  for ($i = 0; $i < $num_rows; $i++) {
    if ($i % 2 == 0) {
      echo "<tr class=\"odd\">";
    } else {
      echo "<tr class=\"even\">";
    }
    $row = $result->fetch_assoc();
    echo '<td><input type="checkbox" name="meta_data[]" value="'.$row[$primary_key].'"></td>';
    $j = 0;
    while ($j < $result->field_count) {
      $meta = $result->fetch_field_direct($j);
      if ($row[$meta->name] == "done"){
        echo '<td headers="' . $meta->name . '" class="metadata"><center><img src="images/done.png" width="20" height="20"></center></td>';
      } else {
        echo '<td headers="' . $meta->name . '" class="metadata"><center>' . $row[$meta->name] . '</center></td>';
      }
      $j++;
    }
    echo "</tr>";
  }
  echo "</table></form>";
}
?>
  
<?php //Authenticate Login
//function atlas_authenticate() {
  //if($_SESSION['user_is_logged_in'] == false){
  //  header('Location: https://geco.iplantcollaborative.org/modupeore17/atlas/index.php');
//  }
//}
?>
<?php //Sequence Page function
function metavw_display($action, $result, $primary_key) {
  $num_rows = $result->num_rows;
  echo '<br><table class="metadata"><tr>';
  echo '<th align="left" width=40pt bgcolor="white"><font size="2" color="red">Select All</font><input type="checkbox" id="selectall" onClick="selectAll(this)" /></th>';
  $meta = $result->fetch_field_direct(0); echo '<th class="metadata" id="' . $meta->name . '">Library id</th>';
  $meta = $result->fetch_field_direct(1); echo '<th class="metadata" id="' . $meta->name . '">Line</th>';
  $meta = $result->fetch_field_direct(2); echo '<th class="metadata" id="' . $meta->name . '">Species</th>';
  $meta = $result->fetch_field_direct(3); echo '<th class="metadata" id="' . $meta->name . '">Tissue</th>';
  $meta = $result->fetch_field_direct(4); echo '<th class="metadata" id="' . $meta->name . '">Total reads</th>';
  $meta = $result->fetch_field_direct(5); echo '<th class="metadata" id="' . $meta->name . '">Mapped reads</th>';
  $meta = $result->fetch_field_direct(6); echo '<th class="metadata" id="' . $meta->name . '">Genes</th>';
  $meta = $result->fetch_field_direct(7); echo '<th class="metadata" id="' . $meta->name . '">Isoforms</th>';
  $meta = $result->fetch_field_direct(8); echo '<th class="metadata" id="' . $meta->name . '">Variants</th>';
  $meta = $result->fetch_field_direct(9); echo '<th class="metadata" id="' . $meta->name . '">SNPs</th>';
  $meta = $result->fetch_field_direct(10); echo '<th class="metadata" id="' . $meta->name . '">INDELs</th>';
  $meta = $result->fetch_field_direct(11); echo '<th class="metadata" id="' . $meta->name . '">Notes</th></tr>';

  for ($i = 0; $i < $num_rows; $i++) {
    if ($i % 2 == 0) {
      echo "<tr class=\"odd\">";
    } else {
      echo "<tr class=\"even\">";
    }
    $row = $result->fetch_assoc();
    echo '<td><input type="checkbox" name="meta_data[]" value="'.$row[$primary_key].'"></td>';
    $j = 0;
    while ($j < $result->field_count) {
      $meta = $result->fetch_field_direct($j);
      echo '<td headers="' . $meta->name . '" class="metadata"><center>' . $row[$meta->name] . '</center></td>';
      $j++;
    }
    echo "</tr>";
  }
  echo "</table></form>";
}
?>

