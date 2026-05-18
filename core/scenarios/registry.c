/* Scenario registry. Each scenario .c file provides a registration entry-point
 * (gr_scenario_register_<name>) that pushes its gr_scenario_t descriptor onto
 * a static table at startup. */

#include "grlite.h"

#include <string.h>

#define GR_MAX_SCENARIOS 32

static const gr_scenario_t* g_scenarios[GR_MAX_SCENARIOS];
static int g_scenario_count = 0;
static int g_initialized = 0;

void gr_scenario_register(const gr_scenario_t* s) {
    if (!s) return;
    if (g_scenario_count >= GR_MAX_SCENARIOS) return;
    g_scenarios[g_scenario_count++] = s;
}

const gr_scenario_t* gr_scenario_find(const char* name) {
    if (!name) return NULL;
    for (int i = 0; i < g_scenario_count; i++) {
        const gr_scenario_t* s = g_scenarios[i];
        if (s && s->name && strcmp(s->name, name) == 0) return s;
    }
    return NULL;
}

/* Each new scenario .c file declares its registration function here and calls
 * it from gr_scenarios_init. */
extern void gr_scenario_register_wave_pulse(void);
extern void gr_scenario_register_static_source(void);

void gr_scenarios_init(void) {
    if (g_initialized) return;
    g_initialized = 1;
    gr_scenario_register_wave_pulse();
    gr_scenario_register_static_source();
}
