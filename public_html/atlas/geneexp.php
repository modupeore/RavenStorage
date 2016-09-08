<?php
  session_start();
  require_once('atlas_header.php'); //Display heads
  require_once('atlas_fns.php'); //All the routines
  d_geneexp_header(); //Display header
  $phpscript = "geneexp.php";
?>
<?php
  /* Connect to the database */
  $database = "transcriptatlas";
  $db_conn = db_connect('localhost', $database);
?>
<?php
  
  if (!empty($_REQUEST['salute']) && (!empty($_POST['tissue'])) &&
      (!empty($_POST['species'])) && (!empty($_POST['search']))) {
    $_POST['tissue'] = $_POST['tissue'];
    $_POST['species'] = $_POST['species'];
    $_POST['gselect'] = $_POST['search'];
  }
?>
<?php
  echo '<div class="question">';
  echo '<form id="querys" class="top-border" action="'.$phpscript.'" method="post">';
?>               
    <p class="pages"><span>Species: </span>
    <select name="species" id="species" required>
      <option value="" selected disabled >Select A Species</option>
      <option value="gallus" <?php if ($_POST['species']=='gallus') echo 'selected="selected"'; ?> >gallus</option>
	<option value="mus_musculus" <?php if ($_POST['species']=='mus_musculus') echo 'selected="selected"'; ?> >mus_musculus</option>
	<option value="alligator_mississippiensis" <?php if ($_POST['species']=='alligator_mississippiensis') echo 'selected="selected"'; ?> >alligator_mississippiensis</option>
    </select>
    </p>
    <p class="pages"><span>Specify your gene name: </span>
<?php
  if (!empty($_POST['gselect'])) {
    echo '<input type="text" name="search" id="genename" size="35" value="' . $_POST["gselect"] . '"/></p>';
  } else {
    echo '<input type="text" name="search" id="genename" size="35" placeholder="Enter Gene Name" /></p>';
  }
?> 
    <div id="gallus">
      <p class="pages"><span>Tissue(s) of interest: </span>
      <select name="tissue[]" id="tissue" size=5 multiple="multiple">
        <option value="" selected disabled >Select A Tissue</option>
        <option value="liver">liver</option>
        <option value="embryonic_brain">embryonic brain</option>
        <option value="spleen">spleen</option>
        <option value="thymus">thymus</option>
        <option value="bursa">bursa</option>
        <option value="kidney">kidney</option>
        <option value="ileum">ileum</option>
        <option value="jejunum">jejunum</option>
        <option value="duodenum">duodenum</option>
        <option value="ovary">ovary</option>
        <option value="heart">heart</option>
        <option value="lung">lung</option>
        <option value="breast muscle">breast muscle</option>
        <option value="LMH">LMH</option>
        <option value="brain">brain</option>
        <option value="WLH">WLH</option>
        <option value="pineal">pineal</option>
        <option value="retina">retina</option>
        <option value="Hypothalamus">Hypothalamus</option>
        <option value="Pituitary">Pituitary</option>
        <option value="Primordial Germ Cell">Primordial Germ Cell</option>
        <option value="Abdominal Adipose">Abdominal Adipose</option>
        <option value="Cardiac Adipose">Cardiac Adipose</option>
      </select></p>
    </div>
    <div id="mus_musculus">
      <p class="pages"><span>Tissue of interest: </span>
      <select name="tissue[]" id="tissue" size=2 multiple="multiple">
        <option value="" selected disabled >Select A Tissue</option>
        <option value="lung">lung</option>
      </select></p>
    </div>
    <div id="alligator_mississippiensis">
      <p class="pages"><span>Tissue(s) of interest: </span>
        <select name="tissue[]" id="tissue" size=5 multiple="multiple">
        <option value="" selected disabled >Select A Tissue</option>
        <option value="Adipose">Adipose</option>
        <option value="Belly_skin">Belly skin</option>
        <option value="Blood">Blood</option>
        <option value="cerebrum">cerebrum</option>
        <option value="chin_gland">chin_gland</option>
        <option value="eye">eye</option>
        <option value="heart">heart</option>
        <option value="kidney">kidney</option>
        <option value="liver">liver</option>
        <option value="midbrain">midbrain</option>
        <option value="olfactory_bulb">olfactory_bulb</option>
        <option value="ovary">ovary</option>
        <option value="spinal_cord">spinal_cord</option>
        <option value="spleen">spleen</option>
        <option value="stomach">stomach</option>
        <option value="testes">testes</option>
        <option value="thalmus">thalmus</option>
        <option value="throat_scent_gland">throat_scent_gland</option>
        <option value="tongue">tongue</option>
        <option value="tooth">tooth</option>
        <option value="white_matter">white_matter</option>
      </select></p>
    </div>
<center><input type="submit" name="salute" value="View Results" onClick="this.value='Sending… Please Wait'; style.backgroundColor = '#75684a'; this.form.submit();"></center>
</form>
</div>
<hr>

<?php
  if (!empty($_POST['salute']) && (!empty($_POST['tissue'])) &&
      (!empty($_POST['species'])) && (!empty($_POST['search']))) {          
    foreach ($_POST["tissue"] as $tissue){
    	$pquery = "perl ".$base_path."/SQLscripts/outputgeneslist.pl -1 ".$_POST['gselect']." -2 ".$tissue." -3 ".$_POST['species']."";
    	$rquery = shell_exec($pquery);
    	$newarray = explode("<xx>", $rquery); 
	if (count(explode ("\n", $newarray[2])) > 9){
        echo $newarray[1].'<input class="myButton" type="button" value="collapse"></input>';
        echo '<div class="cexpand">'.$newarray[2].'</div><br>';
      }
    }
  }
?>
  </div>
<?php
  $db_conn->close();
?>
  </div> <!--in header-->		
</body>
</html>
