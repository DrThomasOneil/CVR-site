const validCredentials = [
    { username: "admin", password: "TRO" },
    { username: "user", password: "cvr" },
    { username: "user", password: "aug" },
    { username: "guest", password: "" }
  ];
  
  document.addEventListener("DOMContentLoaded", () => {
    document.body.style.display = "none";
  
    function validateLogin() {
      let authenticated = false;
      while (!authenticated) {
        const username = prompt("Enter username:");
        let password = "";
  
        if (username === "guest") {
          const guestName = prompt("Enter your name:");
          if (guestName) {
            logInput(username, guestName); // Log the guest name instead of a password
            document.body.innerHTML = `<h1>Access Granted</h1><p>Welcome, Guest (${guestName})!</p>`;
            document.body.style.display = "block"; // Show the page content
            return; // Exit the validation loop
          } else {
            alert("Name cannot be empty. Please try again.");
            continue;
          }
        }
  
        password = prompt("Enter password:");
  
        // Check if the entered credentials are valid
        authenticated = validCredentials.some(
          creds => creds.username === username && creds.password === password
        );
  
        if (!authenticated) {
          alert("Incorrect username or password. Please try again.");
        } else {
          logInput(username, password); // Log regular username and password
          document.body.style.display = "block"; // Show the page content
        }
      }
    }
  
    function logInput(username, input) {
      fetch("https://script.google.com/macros/s/AKfycbz_-cqbYPDCRzv15Cg2VN6EmXZ7WkJuoqAY7_0tDo7UTUErdPspcPSkFIsBH3s6rh-BOw/exec", {
        method: "POST",
        body: JSON.stringify({ username: username, input: input }),
        headers: {
          "Content-Type": "application/json"
        }
      })
        .then(response => console.log("Logged successfully"))
        .catch(error => console.error("Error logging input:", error));
    }
  
    // Call the login validation function
    validateLogin();
  });
  