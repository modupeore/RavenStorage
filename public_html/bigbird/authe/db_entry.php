<?php
  session_start();
  require_once('db_functions.php'); //All the routines
  db_authenticate(); //Authenticate login information
  display_header(); //Display header

  //$phpscript = "db_entry.php"; //apparently it works even when commented out *weird*
?>
<?php
  /* Connect to the database */
  $database = "transcriptatlas";
  $table = "bird_libraries";
  $db_conn = db_connect('localhost', $database);
  if (isset($_POST['accept']))
  {
    $birdid = mysqli_real_escape_string($db_conn, htmlentities($_POST["bird_id"]));
    $line = mysqli_real_escape_string($db_conn, htmlentities($_POST["line"]));
    $tissue = mysqli_real_escape_string($db_conn, htmlentities($_POST["tissue"]));
    $indexs = mysqli_real_escape_string($db_conn, htmlentities($_POST["index_"]));
    $scientist = mysqli_real_escape_string($db_conn, htmlentities($_POST["scientist"]));
    $notes = mysqli_real_escape_string($db_conn, htmlentities($_POST["notes"]));
    $species = $_POST["species"];
    $method = $_POST["method"];
    $chip_result = $_POST["chip_result"];
  }
?>
<?php
/* Rest of It*/
$query = "SHOW KEYS FROM $table WHERE Key_name = 'PRIMARY'";
$result = $db_conn->query($query);
// create array of primary keys
$primary_key = null;
while ($row = $result->fetch_assoc()) {
	$primary_key[] = $row['Column_name'];
}
if (!$primary_key) {
	echo "<div><strong>NO PRIMARY KEY!\nFUNCTIONS WILL NOT WORK AS EXPECTED!</strong></div>";
}
// get total number of rows in table
$query = "SELECT * FROM $table";
$all_rows = $db_conn->query($query);
$total_rows = $all_rows->num_rows;

echo '<form action="'.$phpscript.'" method="post">';
echo '<table><tr><td align="right"><input type="submit" name="logout" value="logout"/></td></tr></table><br><br>';
?>
<table>
<tr><td class="border"></td>
<th class="border" width=6%><strong><font size=2px>bird_id</font></strong></th>
<th class="border"><strong><font size=2px>species</font></strong></th>
<th class="border" width=10%><strong><font size=2px>line</font></strong></th>
<th class="border" width=10%><strong><font size=2px>tissue</font></strong></th>
<th class="border"><strong><font size=2px>method</font></strong></th>
<th class="border" width=6%><strong><font size=2px>index</font></strong></th>
<th class="border"><strong><font size=2px>chip_result</font></strong></th>
<th class="border" width=10%><strong><font size=2px>scientist</font></strong></th>
<th class="border" width=20%><strong><font size=2px>notes</font></strong></th>
</tr><tr>

<!--create accept functionality-->
<td class="border"><input type="submit" name="accept" value="insert"/></td>
<td class="border"><input type="text" class="forms" name="bird_id"<?php if(!empty($db_conn)){echo 'value="'.$birdid.'"';}?>/></td>	<!--bird_id-->
<td class="border"><select name="species">
	<option value="gallus">gallus</option>
	<option value="bos_taurus" <?php if ($_POST['species']=='bos_taurus') echo 'selected="selected"'; ?> >bos_taurus</option>
	<option value="mus_musculus" <?php if ($_POST['species']=='mus_musculus') echo 'selected="selected"'; ?> >mus_musculus</option>
	<option value="meleagris_gallopavo" <?php if ($_POST['species']=='meleagris_gallopavo') echo 'selected="selected"'; ?> >meleagris_gallopavo</option>
	<option value="anatidae" <?php if ($_POST['species']=='anatidae') echo 'selected="selected"'; ?> >anatidae</option>
	<option value="alligator_mississippiensis" <?php if ($_POST['species']=='alligator_mississippiensis') echo 'selected="selected"'; ?> >alligator_mississippiensis</option>
	<option value="homo_sapiens" <?php if ($_POST['species']=='homo_sapiens') echo 'selected="selected"'; ?> >homo_sapiens</option>
	<option value="tardigrade" <?php if ($_POST['species']=='tardigrade') echo 'selected="selected"'; ?> >tardigrade</option>
	</select></td> <!--species-->
<td class="border"><input type="text" class="forms" name="line"<?php if(!empty($db_conn)){echo 'value="'.$line.'"';}?>/></td> <!--line-->
<td class="border"><input type="text" class="forms" name="tissue"<?php if(!empty($db_conn)){echo 'value="'.$tissue.'"';}?>/></td> <!--tissue-->
<td class="border">
        <select name="method">
        <option value="rna_seq">rna_seq</option>
        <option value="micro_seq" <?php if ($_POST['method']=='micro_seq') echo 'selected="selected"'; ?> >micro_seq</option>
        <option value="index_seq" <?php if ($_POST['method']=='index_seq') echo 'selected="selected"'; ?> >index_seq</option>
        <option value="16s" <?php if ($_POST['method']=='16s') echo 'selected="selected"'; ?> >16s</option>
        </select></td> <!--method-->
<td class="border"><input type="text" class="forms" name="index_"<?php if(!empty($db_conn)){echo 'value="'.$indexs.'"';}?>/></td> <!--index-->
<td class="border">
        <select name="chip_result">
        <option value=""></option>
        <option value="accept" <?php if ($_POST['chip_result']=='accept') echo 'selected="selected"'; ?> >accept</option>
        <option value="reject" <?php if ($_POST['chip_result']=='reject') echo 'selected="selected"'; ?> >reject</option>
        </select></td> <!--chip_result-->
<td class="border"><input type="text" class="forms" name="scientist"<?php if(!empty($db_conn)){echo 'value="'.$scientist.'"';}?>/></td> <!--scientist-->
<td class="border"><input type="text" class="forms" name="notes"<?php if(!empty($db_conn)){echo 'value="'.$notes.'"';}?>/></td> <!--note-->		
</tr>
</table>
</form>
<br><hr>
<?php
// delete from table
db_delete($primary_key, $table, $db_conn);
// verify the input
db_accept($phpscript,$table, $db_conn);
//insert into table
db_insert($table, $db_conn);
//
db_reset(); db_logout();
?>
<?php
//create query for DB display
if (!empty($_REQUEST['order'])) {
    // if the sort option was used
    $_SESSION[$table]['sort'] = $_POST['sort'];
    $_SESSION[$table]['dir'] = $_POST['dir'];
    $_SESSION[$table]['num_recs'] = $_POST['num_recs'];

    $terms = explode(",", $_POST['search']);
    $is_term = false;
    foreach ($terms as $term) {
        if (trim($term) != "") {
            $is_term = true;
        }
    }
    $_SESSION[$table]['select'] = $terms;
    $_SESSION[$table]['column'] = $_POST['column'];

    $query = ("SELECT * FROM $table ");
    if ($is_term) {
        $query .= "WHERE ";
    }
    foreach ($_SESSION[$table]['select'] as $term) {
        if (trim($term) == "") {
            continue;
        }
        $query .= $_SESSION[$table]['column'] . " LIKE '%" . trim($term) . "%' OR ";
    }
    $query = rtrim($query, " OR ");
    $query .= " ORDER BY " . $_SESSION[$table]['sort'] . " " . $_SESSION[$table]['dir'];

    $result = $db_conn->query($query);
    $num_total_result = $result->num_rows;
    if ($_SESSION[$table]['num_recs'] != "all") {
        $query .= " limit " . $_SESSION[$table]['num_recs'];
    }
    unset($_SESSION[$table]['txt_query']);
} elseif (!empty($_SESSION[$table]['sort'])) {
    $is_term = false;
    foreach ($_SESSION[$table]['select'] as $term) {
        if (trim($term) != "") {
            $is_term = true;
        }
    }
    $query = ("SELECT * FROM $table ");
    if ($is_term) {
        $query .= "WHERE ";
    }
    foreach ($_SESSION[$table]['select'] as $term) {
        if (trim($term) == "") {
            continue;
        }
        $query .= $_SESSION[$table]['column'] . " LIKE '%" . trim($term) . "%' OR ";
    }
    $query = rtrim($query, " OR ");
    $query .= " ORDER BY " . $_SESSION[$table]['sort'] . " " . $_SESSION[$table]['dir'];
    $result = $db_conn->query($query);
    $num_total_result = $result->num_rows;

    if ($_SESSION[$table]['num_recs'] != "all") {
        $query .= " limit " . $_SESSION[$table]['num_recs'];
    }
} else {
    // if this is the first time, then just order by line and display all rows //default
        $query = "SELECT * FROM $table ORDER BY $primary_key[0] desc limit 10";
}
$result = $db_conn->query($query);
if ($db_conn->errno) {
    echo "<div id=\"index-error\"";
    echo "<span><strong>Error with query.</strong></span>";
    echo "<span><strong>Error number: </strong>$db_conn->errno</span>";
    echo "<span><strong>Error string: </strong>$db_conn->error</span>";
    echo "</div>";
}
$num_results = $result->num_rows;
if (empty($_SESSION[$table]['sort'])) {
    $num_total_result = $num_results;
}
?>
<?php
echo '<form id="sort" class="top-border" action="'.$phpscript.'" method="post">';
?>
<!-- QUERY -->
<div>
    <p class=lspace><span>Search for: </span>
<?php
  if (!empty($_SESSION[$table]['select'])) {
    echo '<input type="text" size="40" name="search" value="' . implode(",", $_SESSION[$table]["select"]) . '"\"/>';
  } else {
    echo '<input type="text" size="40" name="search" placeholder="Enter variable(s) separated by commas (,)"/>';
  }
   
?>
    <span> in </span>
    <select name="column">
      <option value="library_id">library_id</option>
      <?php
        if (empty($_SESSION[$table]['column'])) {
          $_SESSION[$table]['library_id'] = "library_id";
        }
        if ($_SESSION[$table]['column'] == "bird_id") {
          echo '<option selected value="bird_id">bird_id</option>';
        } else {
          echo '<option value="bird_id">bird_id</option>';
        }
        if ($_SESSION[$table]['column'] == "species") {
          echo '<option selected value="species">species</option>';
        } else {
          echo '<option value="species">species</option>';
        }
        if ($_SESSION[$table]['column'] == "line") {
          echo '<option selected value="line">line</option>';
        } else {
          echo '<option value="line">line</option>';
        }
        if ($_SESSION[$table]['column'] == "tissue") {
          echo '<option selected value="tissue">tissue</option>';
        } else {
          echo '<option value="tissue">tissue</option>';
        }
        if ($_SESSION[$table]['column'] == "method") {
          echo '<option selected value="method">method</option>';
        } else {
          echo '<option value="method">method</option>';
        }
        if ($_SESSION[$table]['column'] == "chip_result") {
          echo '<option selected value="chip_result">chip_result</option>';
        } else {
          echo '<option value="chip_result">chip_result</option>';
        }
        if ($_SESSION[$table]['column'] == "scientist") {
          echo '<option selected value="scientist">scientist</option>';
        } else {
          echo '<option value="scientist">scientist</option>';
        }
        if ($_SESSION[$table]['column'] == "notes") {
          echo '<option selected value="notes">notes</option>';
        } else {
          echo '<option value="notes">notes</option>';
        }
      ?> 
    </select></p>
    <p class=lspace ><span>Sort by:</span>
    <select name="sort">
      <option value="library_id">library_id</option>
      <?php
        if (empty($_SESSION[$table]['sort'])) {
          $_SESSION[$table]['library_id'] = "library_id";
        }
        if ($_SESSION[$table]['sort'] == "bird_id") {
          echo '<option selected value="bird_id">bird_id</option>';
        } else {
          echo '<option value="bird_id">bird_id</option>';
        }
        if ($_SESSION[$table]['sort'] == "species") {
          echo '<option selected value="species">species</option>';
        } else {
          echo '<option value="species">species</option>';
        }
        if ($_SESSION[$table]['sort'] == "line") {
          echo '<option selected value="line">line</option>';
        } else {
          echo '<option value="line">line</option>';
        }
        if ($_SESSION[$table]['sort'] == "tissue") {
          echo '<option selected value="tissue">tissue</option>';
        } else {
          echo '<option value="tissue">tissue</option>';
        }
        if ($_SESSION[$table]['sort'] == "method") {
          echo '<option selected value="method">method</option>';
        } else {
          echo '<option value="method">method</option>';
        }
        if ($_SESSION[$table]['sort'] == "index_") {
          echo '<option selected value="index_">index</option>';
        } else {
          echo '<option value="index_">index</option>';
        }
        if ($_SESSION[$table]['sort'] == "chip_result") {
          echo '<option selected value="chip_result">chip_result</option>';
        } else {
          echo '<option value="chip_result">chip_result</option>';
        }
        if ($_SESSION[$table]['sort'] == "scientist") {
          echo '<option selected value="scientist">scientist</option>';
        } else {
          echo '<option value="scientist">scientist</option>';
        }
        if ($_SESSION[$table]['sort'] == "date") {
          echo '<option selected value="date">date</option>';
        } else {
          echo '<option value="date">date</option>';
        }
        if ($_SESSION[$table]['sort'] == "notes") {
          echo '<option selected value="notes">notes</option>';
        } else {
          echo '<option value="notes">notes</option>';
        }
      ?> 
    </select> <!if ascending or descending...>
    <select name="dir">
      <option value="asc">ascending</option>
      <?php
        if (empty($_SESSION[$table]['dir'])) {
          $_SESSION[$table]['asc'] = "asc";
        }
        if ($_SESSION[$table]['dir'] == "desc") {
          echo '<option selected value="desc">descending</option>';
        } else {
          echo '<option value="desc">descending</option>';
        }
      ?>
    </select>
    <span>and show</span>
    <select name="num_recs">
      <option value="10">10</option>
      <?php
        if (empty($_SESSION[$table]['num_recs'])) {
          $_SESSION[$table]['num_recs'] = "10";
        }
        if ($_SESSION[$table]['num_recs'] == "20") {
          echo '<option selected value="20">20</option>';
        } else {
          echo '<option value="20">20</option>';
        }
        if ($_SESSION[$table]['num_recs'] == "50") {
          echo '<option selected value="50">50</option>';
        } else {
          echo '<option value="50">50</option>';
        }
        if ($_SESSION[$table]['num_recs'] == "all") {
          echo '<option selected value="all">all</option>';
        } else {
          echo '<option value="all">all</option>';
        }
      ?> 
    </select>
    <span>records.</span>
    <input type="submit" name="order" value="Go"/></p>
  </div>
</form>
<br/>
<div id="db-display" class="top-border">
<?php
echo "<span>" . $num_results . " out of " . $num_total_result . " search results displayed. (" . $total_rows . " total rows)</span>";
db_display($phpscript, $result, $primary_key);
?>
</div>
<?php
$result->free();
$db_conn->close();
?>
</div>		
<?php
echo "</body>";
echo "</html>";
?>

