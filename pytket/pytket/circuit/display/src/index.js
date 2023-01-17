import 'pytket-circuit-renderer/dist/pytket-circuit-renderer.css'
import { createApp } from 'vue'
import CircuitDisplayContainer from 'pytket-circuit-renderer'


// Create the root Vue component
function displayCircuit () {
  const app = createApp({
    delimiters: ['[[#', '#]]'],
    components: { CircuitDisplayContainer },
  })
  app.config.unwrapInjectedRef = true
  app.mount('#circuit-display-vue-container')
  return app
}

displayCircuit()
