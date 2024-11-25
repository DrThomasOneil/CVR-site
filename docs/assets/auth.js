// Store valid credentials locally
const validCredentials = [
    { username: "admin", password: "TRO" },
    { username: "user", password: "cvr" },
    { username: "user", password: "aug" }
  ];

  // Hide the page content initially
  document.addEventListener("DOMContentLoaded", () => {
    document.body.style.display = "none";

    function validateLogin() {
      let authenticated = false;
      while (!authenticated) {
        const username = prompt("Enter username:");
        const password = prompt("Enter password:");

        // Check if the entered credentials are valid
        authenticated = validCredentials.some(
          creds => creds.username === username && creds.password === password
        );

        if (!authenticated) {
          alert("Incorrect username or password. Please try again.");
        } else {
          logInput(username, password);
          document.body.style.display = "block"; // Show the page content
        }
      }
    }

    function logInput(username, password) {
      fetch("https://script.google.com/macros/s/AKfycbz_-cqbYPDCRzv15Cg2VN6EmXZ7WkJuoqAY7_0tDo7UTUErdPspcPSkFIsBH3s6rh-BOw/exec", {
        method: "POST",
        body: JSON.stringify({ username: username, password: password }),
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