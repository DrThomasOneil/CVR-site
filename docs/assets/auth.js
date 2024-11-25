fetch('https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/docs/assets/cred.js')
    .then(response => response.json())
    .then(data => {
        const validCredentials = data;

        const username = prompt("Enter username:");
        const password = prompt("Enter password:");

        const isValid = validCredentials.some(
            creds => creds.username === username && creds.password === password
        );

        if (isValid) {
            logInput(username, password);
            document.body.innerHTML = "<h1>Access Granted</h1><p>Welcome, " + username + "!</p>";
        } else {
            document.body.innerHTML = "Access Denied. Incorrect username or password.";
        }
    })
    .catch(error => {
        console.error("Error loading credentials:", error);
        document.body.innerHTML = "Error loading authentication data.";
    });

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