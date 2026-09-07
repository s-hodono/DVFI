/*
   TTL modulator for ARFI (Interrupt-based)
   Author: Shota Hodono (modified)
   Date: Oct 2025
*/

const int outputPin = 3;  // TTL output signal pin
const int inputPin  = 2;  // TTL input signal pin (must be interrupt-capable)
const int ledPin    = 4;  // LED pin

const int nSlice = 1;     // number of slices/volume TR

unsigned long delayBeforeOutput = 6;   // 2ms
unsigned long outputDuration    = 1;   // ms
unsigned long maxNoTTLTime      = 50000; // ms

volatile bool ttlReceived = false;
volatile unsigned long lastSignalTime = 0;
bool processTTL = true;

int ttlCounter = 0;
int cycleCounter = 0;

void setup() {
  Serial.begin(9600);
  while (!Serial);

  pinMode(inputPin, INPUT);
  pinMode(outputPin, OUTPUT);
  pinMode(ledPin, OUTPUT);
  digitalWrite(outputPin, LOW);
  digitalWrite(ledPin, LOW);

  attachInterrupt(digitalPinToInterrupt(inputPin), onTTL, RISING);
  Serial.println("Interrupt mode initialized.");
  Serial.print("number of slices:                                                                 ");
  Serial.println(nSlice);
}

void loop() {
  if (ttlReceived) {
    ttlReceived = false;  // reset flag

    ttlCounter++;

    if (processTTL) {
      delay(delayBeforeOutput);
      digitalWrite(outputPin, HIGH);
      delay(outputDuration);
      digitalWrite(outputPin, LOW);
      Serial.println(ttlCounter);
    }

    // Handle slice grouping
    if (ttlCounter == nSlice) {
      ttlCounter = 0;
      cycleCounter++;

      if (cycleCounter == 1) processTTL = false;
      else if (cycleCounter == 2) {
        processTTL = true;
        cycleCounter = 0;
      }
    }
  }

  // Timeout watchdog
//  if (millis() - lastSignalTime > maxNoTTLTime) {
//    Serial.println("No TTL signal, system stopping...");
//    while (true);
//  }
}

void onTTL() {
  unsigned long now = millis();
//  if (now - lastSignalTime > 100) {  // debounce (ms)
    lastSignalTime = now;
    ttlReceived = true;
//  }
}
