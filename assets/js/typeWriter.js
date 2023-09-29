var i = 0;
var txt = 'Passionate about changing the world with technology and exploring interdisciplinary boundaries.';
var speed = 50;

function typeWriter() {
  if (i < txt.length) {
    document.getElementById("subtitle").innerHTML += txt.charAt(i);
    i++;
    setTimeout(typeWriter, speed);
  }
}

window.onload = typeWriter;
