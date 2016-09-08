<?php
//*** DISPLAY FUNCTION ***/
require_once('db_functions.php');
?>
<?php
/*LOGIN AND LOGOUT INFO */
function db_authenticate() {
    if($_SESSION['user_is_logged_in'] == false){
	header('Location: https://geco.iplantcollaborative.org/modupeore17/bigbird/index.php');
    }
}

function db_logout() {
    if (!empty($_REQUEST['logout'])) {
      session_unset();
      session_destroy();
      header('Location: https://geco.iplantcollaborative.org/modupeore17/bigbird/index.php');
    }
}

function db_reset() {
    if (!empty($_REQUEST['reset'])) {
      header('Location: https://geco.iplantcollaborative.org/modupeore17/bigbird/db_entry.php');
    }
}
?>
<?php
function display_majorheader() {
  /*for the index.php page */
    echo "<head>";
    echo "<style>";
    echo 'h1 {
    display:block;
    margin-top:2em;
    font-size:28pt; 
    text-align:center;
    font-family: "calibri", arial,sans-serif;
    font-weight:bold;
    }';
    echo '</style>';
    echo "<meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    echo "<title>Bigbird_DB-Login</title>";
    echo '<img src="title.png" align="middle" alt="BigBird Libraries" style="width:100%;height:100px">';
    echo "</head><body>";
    echo "<div id=\"container\"><br><br>";
}
?>
<?php
function display_header() {
    echo "<head>";
    echo "<meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    echo "<title>Bigbird_DB-Entry</title>";
    echo '<link rel="stylesheet" type="text/css" href="stylebird.css">';
    echo '<img src="title.png" align="middle" alt="BigBird Libraries" style="width:100%;height:100px;">';
    echo "</head><body>";
    echo "<div id=\"container\">";   
}
?>
<?php
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
<?php
//*****DATABASE fUNCTIONALITY*****//
function db_delete($primary_key, $table, $db_conn) {
    if (!empty($_REQUEST['delete'])) {
        foreach ($_POST['check'] as $entry) {
            $entry_arr = explode(",", $entry);
            $query = "DELETE FROM $table WHERE ";
            $i = 1;
            foreach ($primary_key as $pk) {
                $query = $query . "$pk='$entry_arr[$i]' AND ";
                $i++;
            }
            $query = rtrim($query, " AND ");
            $result = $db_conn->query($query);
        }
    }
}
?>
<?php
function db_accept($phpscript, $table, $db_conn) {
    if (!empty($_REQUEST['accept'])) {
        echo '<form action="'.$phpscript.'" method="post">';
  ?>
	  <table>
	    <tr>
		<td align="right" class="lines"><input type="submit" name="reset" value="reject"/>
		  <input type="submit" name="verified" value="accept"/>
		</td>
		<td>
		  <table border=1>
		    <tr>
			<th class="lines"><strong><font size=2px>bird_id</font></strong></th>
			<th class="lines"><strong><font size=2px>species</font></strong></th>
			<th class="lines"><strong><font size=2px>line</font></strong></th>
			<th class="lines"><strong><font size=2px>tissue</font></strong></th>
			<th class="lines"><strong><font size=2px>method</font></strong></th>
			<th class="lines"><strong><font size=2px>index</font></strong></th>
			<th class="lines"><strong><font size=2px>chip_result</font></strong></th>
			<th class="lines"><strong><font size=2px>scientist</font></strong></th>
			<th class="lines"><strong><font size=2px>notes</font></strong></th>
		    </tr>
		    <tr>
		  <?PHP
		    echo '<td class="lines"><input type="hidden" name="bird_id" value="'.$_POST["bird_id"].'"/>'.$_POST["bird_id"].'</td>';	// bird_id
		    echo '<td class="lines"><input type="hidden" name="species" value="'.$_POST["species"].'"/>'.$_POST["species"].'</td>';	// species
		    echo '<td class="lines"><input type="hidden" name="line" value="'.$_POST["line"].'"/>'.$_POST["line"].'</td>';	// line
		    echo '<td class="lines"><input type="hidden" name="tissue" value="'.$_POST["tissue"].'"/>'.$_POST["tissue"].'</td>';	// tissue
		    echo '<td class="lines"><input type="hidden" name="method" value="'.$_POST["method"].'"/>'.$_POST["method"].'</td>';	// method
		    echo '<td class="lines"><input type="hidden" name="index_" value="'.$_POST['index_'].'"/>'.$_POST['index_'].'</td>';	// index
		    echo '<td class="lines"><input type="hidden" name="chip_result" value="'.$_POST['chip_result'].'"/>'.$_POST['chip_result'].'</td>';	// bird_id
		    echo '<td class="lines"><input type="hidden" name="scientist" value="'.$_POST['scientist'].'"/>'.$_POST['scientist'].'</td>';	// bird_id
		    echo '<td class="lines"><input type="hidden" name="notes" value="'.$_POST['notes'].'"/>'.$_POST['notes'].'</td>';	// bird_id
		    echo '</tr></table></td></tr></table></form>';
    }
}
?>
		    
<?php
function db_insert($table, $db_conn) {
    if (!empty($_REQUEST['verified'])) {
        $query = "SELECT * FROM $table";
        $result = $db_conn->query($query);
        $query = "INSERT INTO $table VALUES (";
        $i = 0;
        while ($i < $result->field_count) {
            $meta = $result->fetch_field_direct($i);
            if ($meta->flags & 512 || $meta->flags & 1024) {
                // if auto_increment or timestamp, insert as null
                $query = $query . "NULL" . ",";
            } elseif ($meta -> flags & 2) {
                $db_result = $db_conn->query("select max(library_id) max from bird_libraries");
		    if ($db_result->num_rows > 0) {while($row = $db_result->fetch_assoc()) {$maxnumber = $row["max"];}}
		    $maxnumber = $maxnumber+1;
                $query = $query . "'".$maxnumber."',";
		} elseif (($meta->type == 253 || $meta->type == 254 || $meta->type == 10
              || $meta->type == 11) && !empty($_POST[$meta->name])) {
                // if it's a text field or date, add single quotes
                $query = $query . "'" . $_POST[$meta->name] . "'" . ",";
            } elseif ($meta->type == 12){
		    $query = $query. "'". date('Y-m-d H:i:s')."',";
		}
		elseif (empty($_POST[$meta->name])) {
                $query = $query . "NULL" . ",";
            } else {
                $query = $query . $_POST[$meta->name] . ",";
            }
            $i++;
        }
        $query = rtrim($query, ",");
        $query = $query . ")";

        $result = $db_conn->query($query);
        echo '<div id="insert-status"';
        if (!$result) {
            echo "<span><strong>Insert unsuccessful.</strong></span>";
            echo "<span><strong>Query: </strong>$query</span>";
            echo "<span><strong>Errormessage: </strong>" . $db_conn->error . "</span>";
        } else {
            echo "<span><strong>Insert successful.</strong></span>";
		echo '<span><strong>Library ID is '.$maxnumber.'.</strong></span>';
		echo "<span>" . $db_conn->affected_rows . " row inserted into $table </span>";
        }
        echo '</div>';
    }
}
?>
<?php
function db_display($action, $result, $primary_key) {
    $num_rows = $result->num_rows;
    echo '<form action="' . $action . '" method="post">';
    echo '<table class="result"><tr><td></td>';
    $meta = $result->fetch_field_direct(0); echo '<th class="result" id="' . $meta->name . '">' . library_id . '</th>';
    $meta = $result->fetch_field_direct(1); echo '<th class="result" id="' . $meta->name . '">' . bird_id . '</th>';
    $meta = $result->fetch_field_direct(2); echo '<th class="result" id="' . $meta->name . '">' . species . '</th>';
    $meta = $result->fetch_field_direct(3); echo '<th class="result" id="' . $meta->name . '">' . line . '</th>';
    $meta = $result->fetch_field_direct(4); echo '<th class="result" id="' . $meta->name . '">' . tissue . '</th>';
    $meta = $result->fetch_field_direct(5); echo '<th class="result" id="' . $meta->name . '">' . method . '</th>';
    $meta = $result->fetch_field_direct(6); echo '<th class="result" id="' . $meta->name . '">' . index . '</th>';
    $meta = $result->fetch_field_direct(7); echo '<th class="result" id="' . $meta->name . '">' . chip_result . '</th>';
    $meta = $result->fetch_field_direct(8); echo '<th class="result" id="' . $meta->name . '">' . scientist . '</th>';
    $meta = $result->fetch_field_direct(9); echo '<th class="result" id="' . $meta->name . '">' . date . '</th>';
    $meta = $result->fetch_field_direct(10); echo '<th class="result" width=20% id="' . $meta->name . '">' . notes . '</th>';
    echo '</tr>';
    for ($i = 0; $i < $num_rows; $i++) {
	  if ($i % 2 == 0) {
            echo "<tr class=\"even\">";
        } else {
            echo "<tr class=\"odd\">";
        }
	  $row = $result->fetch_assoc();
        $pk_values_string = "";
        foreach ($primary_key as $pk) {
            $pk_values_string = "$pk_values_string,$row[$pk]";
        }
        $pk_values_string = ltrim($pk_values_string, ",");
        echo "<td class=\"border lock\" id=\"pk_$pk_values_string\"></td>";
        $j = 0;
        while ($j < $result->field_count) {
            $meta = $result->fetch_field_direct($j);
            if (in_array($meta->name, $primary_key)) {
                echo '<td headers="' . $meta->name . '" class="border primary_key edit_off pk_' . $pk_values_string . '"><center>' . $row[$meta->name] . '</center></td>';
            } else {
                echo '<td headers="' . $meta->name . '" class="border edit_off pk_' . $pk_values_string . '"><center>' . $row[$meta->name] . '</center></td>';
            }
            $j++;
        }
        echo "</tr>";
    }
    echo '</table></form>';
}
?>

