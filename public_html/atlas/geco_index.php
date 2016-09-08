<?php
 require_once('atlas_fns.php');
 $phpscript = "index.php";
?>
<!DOCTYPE html>
  <head>
    <link rel="STYLESHEET" type="text/css" href="stylefile.css">
    <link rel="STYLESHEET" type="text/css" href="mainPage.css">
    <link rel="icon" type="image/ico" href="images/atlas_ico.png"/>
    <title>
      TranscriptAtlas
    </title>
  </head>
  <body>
  <div class="allofit">
    <table>
      <tr>
        <td width=30px></td>
        <td align="center">
          <a href="menu.php"><img src="images/atlas_main.png" alt="Transcriptome Atlas" ></a>
        </td>
      </tr>
    </table>
    <center>
      <div class="container">
        <div class="contact">
<?php
        echo '<form method="post" action="' . $phpscript . '" name="loginform">';
        echo '<table width="80%">';
  class autoLogin 
{
    /**
     * Admin Name
     */
    public $admin = "";

    /**
     * @var bool Login status of user
     */
    private $user_is_logged_in = false;

    /**
     * @var string System messages, likes errors, notices, etc.
     */
    public $feedback = "";

    /**
     * Subdirectories Path
     */
    public $subdirectories = "";
	
	/**
	 * OAuth2 Secret and Key
	 */
	private $secret = '';
	private $key = '';
	
	public function __construct()
	{
	/**
     	  * Configuration file loading
     	 */ 
	$ini = parse_ini_file('/opt/apache2/frankdec/config.ini.php', true);
    	$this->admin = $ini['login']['admin'];
    	$this->subdirectories = $ini['filepaths']['subdirectories'];
	$this->secret = $ini['login']['secret'];
	$this->key = $ini['login']['key'];

        if ($this->performMinimumRequirementsCheck()) {
            $this->runApplication();
        }
    }

    /**
     * Performs a check for minimum requirements to run this application.
     * Does not run the further application when PHP version is lower than 5.3.7
     * Does include the PHP password compatibility library when PHP version lower than 5.5.0
     * (this library adds the PHP 5.5 password hashing functions to older versions of PHP)
     * @return bool Success status of minimum requirements check, default is false
     */
    private function performMinimumRequirementsCheck()
    {
        if (version_compare(PHP_VERSION, '5.3.7', '<')) {
            echo "Sorry, Simple PHP Login does not run on a PHP version older than 5.3.7 !";
        } elseif (version_compare(PHP_VERSION, '5.5.0', '<')) {
            require_once("libraries/password_compatibility_library.php");
            return true;
        } elseif (version_compare(PHP_VERSION, '5.5.0', '>=')) {
            return true;
        }
        // default return
        return false;
    }	 #PHP version
	public function runApplication(){
		// start the session, always needed!
		session_start();
		// check for possible user interactions (login with session/post data or logout)
		$this->performUserLoginAction();
		// show "page", according to user's login status
		if ($this->getUserLoginStatus()) {
			$this->showPageLoggedIn();
			echo '<input type="hidden" name="authenticate" value="true"/>';
		} else {
			echo '<input type="hidden" name="authenticate" value="false"/>';
			$this->showPageLoginForm();
		}
	}
	public function getUserLoginStatus()
	{
        return $this->user_is_logged_in;
	}
	private function performUserLoginAction()
	{
        if (isset($_GET["action"]) && $_GET["action"] == "logout") {
            //$this->doLogout();
	    echo "";
        } elseif (!empty($_SESSION['user_name']) && ($_SESSION['user_is_logged_in'])) {
			$this->doLoginWithSessionData();
        } elseif (isset($_POST["login"])) {
            $this->doLoginWithPostData();
        }
	}
	private function doLoginWithSessionData()
	{
        $this->user_is_logged_in = true; 
	$_SESSION['user_is_logged_in'] = true;
	}

    /**
     * Process flow of login with POST data
     */
    private function doLoginWithPostData() {
		if ($this->checkLoginFormDataNotEmpty()) {
		
			/**
			 * Generate access token using Agave API
			 */
			 
			$ch = curl_init();

			$pf = "grant_type=password&username=".$_POST['user_name']."&password=".$_POST['user_password']."&scope=PRODUCTION";
			$key_and_secret = $this->key.":".$this->secret;
			$encoding = "Content-Type: application/x-www-form-urlencoded";
			
			curl_setopt($ch, CURLOPT_URL, "https://agave.iplantc.org/token");
			curl_setopt($ch, CURLOPT_SSL_VERIFYPEER, 0);
			curl_setopt($ch, CURLOPT_CUSTOMREQUEST, "POST");
			curl_setopt($ch, CURLOPT_POSTFIELDS, $pf);
			curl_setopt($ch, CURLOPT_USERPWD, $key_and_secret);
			curl_setopt($ch, CURLOPT_ENCODING, $encoding);
			curl_setopt($ch, CURLOPT_RETURNTRANSFER, 1);
			
			$response = curl_exec($ch);
			$response_arr = json_decode($response, true);
			$_SESSION['access_token'] = $response_arr['access_token'];
			
			curl_close ($ch);
			/**
			 * Get Login info using access token
			 */
			 
			$ch = curl_init();
			$data = array("Authorization: Bearer ".$_SESSION['access_token']);
			$url = "https://agave.iplantc.org:443/profiles/v2/".$_POST['user_name'];
			curl_setopt($ch, CURLOPT_URL, $url);
			curl_setopt($ch, CURLOPT_HTTPHEADER, $data);
			curl_setopt($ch, CURLOPT_RETURNTRANSFER, 1);
			curl_setopt($ch, CURLOPT_VERBOSE, 1);
						
			$response = curl_exec($ch);
			curl_close ($ch);
			$response_arr = json_decode($response, true);
			
			if($response_arr['status'] == "success"){
				$user_info = $response_arr['result'];
				$_SESSION['user_name'] = $user_info['username'];
				$_SESSION['user_email'] = $user_info['email'];
				$this->user_is_logged_in = true;
				$_SESSION['user_is_logged_in'] = true;
			} else {
				$this->feedback = "Invalid Username or Password ";//.$_SESSION['access_token'];
			}
			
			
		}
    }
	private function checkLoginFormDataNotEmpty(){
        if (!empty($_POST['user_name']) && !empty($_POST['user_password'])) {
            return true;
        } elseif (empty($_POST['user_name'])) {
            $this->feedback = "Username field was empty.";
        } elseif (empty($_POST['user_password'])) {
            $this->feedback = "Password field was empty.";
        }
        // default return
        return false;
    }
	private function showPageLoggedIn(){
		if (isset($_SESSION['LAST_ACTIVITY']) && (time() - $_SESSION['LAST_ACTIVITY'] > 1800)) {
    			// last request was more than 30 minutes ago
    			session_unset();     // unset $_SESSION variable for the run-time 
    			session_destroy();   // destroy session data in storage
		}

		$_SESSION['LAST_ACTIVITY'] = time(); // update last activity time stamp
		$_SESSION['SESSION_TIMEOUT'] = time(); //set timeout for session data

		header('Location: menu.php');
	}
	private function showPageLoginForm(){
		echo "<tr><th colspan=100%><center><h3>TranscriptAtlas Portal Login!</h3></center></th></tr>";
		echo '<tr><td width="160">';
		echo '<label for="login_input_username">Username (or email) : </label></td> ';
		echo '<td><input id="login_input_username" class="contact" type="text" name="user_name" required /> </td></tr>';
		echo '<tr><td><label for="login_input_password">Password : </label> </td>';
		echo '<td><input id="login_input_password" class="contact" type="password" name="user_password" required /></td></tr>';
		echo '<tr><td colspan="100%"><center><input type="submit"  class="contact" name="login" value="Log in" /><center>';
		echo '</td></tr></table>
                  </form></div></div>';
            if ($this->feedback) {
              echo "<font color=red>".$this->feedback."</font><br/><br/>";
            }
            echo "<br><br>";
	}
}

// run the application
$application = new autoLogin();

?>
