<?php
  session_start();
  require_once('atlas_header.php'); //Display heads
  require_once('atlas_fns.php'); //All the routines
  d_kaks_header();
  $phpscript = "2kaks.php";
  $json = file_get_contents('kaks/rossillinoisBM.json');
  $yummy = json_decode($json,true);
  unset($_SESSION['kaksdata']);
?>
  <div class="explain"><p><center>This provides a list of compiled kaks from a subset of our Ross and Illinois libraries.
  <br>
  <br>
  View kaks ratios based on a specific gene of interest.</center></p></div>
<?php

  echo '<form id="query" class="top-border" action="'.$phpscript.'" method="post">'; 
?>

    <p class="pages"><span>Specify your gene name: </span>
      <?php

        if (!empty($_POST['reveal']) && !empty($_POST['select'])) { /* GeneName */
           $_SESSION['select'] = $_POST['select'];
        }
  
          if (!empty($_SESSION['select'])) {
            echo '<input type="text" name="select" id="genename" value="' . $_SESSION["select"] . '"/>';
          } else {
            echo '<input type="text" name="select" id="genename" placeholder="Enter Gene Name" />';
          }
      ?>
      </p><br>
    <center><input type="submit" name="reveal" value="View Results"></center>
  </form></div> <hr>
     <?php
  if (isset($_POST['reveal']) && !empty($_SESSION['select'])) {
    echo '<table><tr><td align="right">';
    foreach ($yummy[$_SESSION['select']] as $sabc) {
      echo '<table class="gened">
            <tr>
              <th class="geneds">Line</th>
              <th class="geneds">Htseq</th>
              <th class="geneds">Cufflinks</th>
            </tr>';
      echo '<tr>
              <td class="gened"><b>Ross</b></td>
              <td class="gened">'.$sabc['Rhtseq'].'</td>
              <td class="gened">'.$sabc['Rcuff'].'</td>
            </tr>';
      echo '<tr>
              <td class="gened"><b>Illinois</b></td>
              <td class="gened">'.$sabc['Ihtseq'].'</td>
              <td class="gened">'.$sabc['Icuff'].'</td>
            </tr>';
      echo '</table><br>';
      break;
    }
    echo '</td></tr><tr></tr><tr><td>';
    echo '<table class="gened">';
    echo '<tr>
            <th class="gened" colspan=5>Chromosomal Location</th>
            <th></th>
            <th class="gened" colspan=5>Ross</th>
            <th></th>
            <th class="gened" colspan=5>Illinois</th>
          </tr>';      
    echo '<tr>
            <th class="geneds" colspan=2>Chr</th>
            <th class="geneds">Start</th>
            <th class="geneds">Stop</th>
            <th class="geneds">Orn</th>
            <th></th>
            <th class="geneds">MLWL</th>
            <th class="geneds">GY</th>
            <th class="geneds">YN</th>
            <th class="geneds">MYN</th>
            <th class="geneds">MA</th>
            <th></th>
            <th class="geneds">MLWL</th>
            <th class="geneds">GY</th>
            <th class="geneds">YN</th>
            <th class="geneds">MYN</th>
            <th class="geneds">MA</th>
          </tr>';
    $number = 0; $list = array(); $line = array();
    foreach ($yummy[$_SESSION['select']] as $sabc) {
      $number += 1;
      
      print '<tr><td class="gened">'.$number.'</td>
            <td class="gened">'.$sabc['chr'].'</td>
            <td class="gened">'.$sabc['start'].
            '</td><td class="gened">'.$sabc['stop'].'</td>
            <td class="gened">'. $sabc['orn'].'</td>
            <td></td>
            <td class="gened">'.$sabc['Rmlwl'].'.</td>
            <td class="gened">'.$sabc['Rgy'].'</td>
            <td class="gened">'.$sabc['Ryn'].'</td>
            <td class="gened">'.$sabc['Rmyn'].'</td>
            <td class="gened">'.$sabc['Rma'].'</td>
            <td></td>
            <td class="gened">'.$sabc['Imlwl'].'</td>
            <td class="gened">'. $sabc['Igy'].'</td>
            <td class="gened">'.$sabc['Iyn'].'</td>
            <td class="gened">'.$sabc['Imyn'].'</td>
            <td class="gened">'.$sabc['Ima']."</td>
          </tr>";

      #Graph image options
      if ($sabc['Rmlwl'] == 50) {$sabc['Rmlwl'] = 0;} if ($sabc['Imlwl'] == 50) {$sabc['Imlwl'] = 0;}
      if ($sabc['Rgy'] == 50) {$sabc['Rgy'] = 0;} if ($sabc['Igy'] == 50) {$sabc['Igy'] = 0;}
      if ($sabc['Ryn'] == 50) {$sabc['Ryn'] = 0;} if ($sabc['Iyn'] == 50) {$sabc['Iyn'] = 0;}
      if ($sabc['Rmyn'] == 50) {$sabc['Rmyn'] = 0;} if ($sabc['Imyn'] == 50) {$sabc['Imyn'] = 0;}
      if ($sabc['Rma'] == 50) {$sabc['Rma'] = 0;} if ($sabc['Ima'] == 50) {$sabc['Ima'] = 0;}
      
      $Rsum = 0; $Rsum = (($sabc['Rmlwl']+$sabc['Rgy']+$sabc['Ryn']+$sabc['Rmyn']+$sabc['Rma'])/5);
      $Isum = 0; $Isum = (($sabc['Imlwl']+$sabc['Igy']+$sabc['Iyn']+$sabc['Imyn']+$sabc['Ima'])/5);
      
      $naming = $_SESSION['select'];
      #pushing to the big array
      array_push($list,$Rsum,$Isum);
    }
    echo "</table>";
    echo '</td></tr></table>';
    #formating the array to be php plot friendly
    $result = count($list); $x = 1; 
    $newlist = $list[0];
    for ($x = 1; $x < $result; $x++) {$newlist .= "," .$list[$x];}

#    exec("Rscript my_rscript.R $newlist "); #$number $_SESSION['select']");
    echo "Rscript my_rscript.R $newlist $number $naming";# $line $_SESSION['select']");


    #exec("Rscript my_rscript.R $_SESSION['kaksdata']");
    #return image
    echo ("<img src='OUTPUT/temp.png />");


    echo '<form id="query" class="top-border" target="_blank" action="bargraph.php" method="post">';
    echo '<center><input type="submit" name="result" value="view bar graph"/></center>';
    echo '</form>';
  }
?>
