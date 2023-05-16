// Script to initialise the circuit renderer app

const { createApp } = Vue;
const circuitDisplayContainer = window["pytket-circuit-renderer"].default;
// Init variables to be shared between circuit display instances
if (typeof window.pytketCircuitDisplays === "undefined") {
    window.pytketCircuitDisplays = {};
}
// Create the root Vue component
const app = createApp({
    delimiters: ['[[#', '#]]'],
    components: { circuitDisplayContainer },
    data () {
      return {
        initRenderOptions: displayOptions,
      }
    }
})
app.config.unwrapInjectedRef = true;
app.mount("#circuit-display-vue-container-"+circuitRendererUid);
window.pytketCircuitDisplays[circuitRendererUid] = app;