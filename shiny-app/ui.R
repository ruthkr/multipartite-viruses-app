# Libraries ------------------------------------------------

# UI Libs
library(shinydashboard)
library(dygraphs)
library(shinyjs)

# UI -------------------------------------------------------

header <- dashboardHeader(
	title = span(tagList(icon("flask"), "Multipartite Viruses"))
)

plot.checkbox <- checkboxGroupInput(
	"plots_checkbox", "Plots",
	choices = c(
		"Time Series" = "time_series",
		"Phase Portrait" = "phase_port",
		"Bifurcation Diagram" = "bifur_diagram",
		"Variation of Gamma" = "gamma_var"
	)
)

plot.inputs <- box(width = 12,
	sliderInput("slider_r1", label = "R1", min = 0, max = 1, step = 0.1, value = 0.1)
)

plot.params <- box(width = 12,
	sliderInput("slider_kappa",  label = "Kappa", min = 0, max = 1, step = 0.1, value = 1),
	sliderInput("slider_gamma1", label = "Gamma", min = 0, max = 1, step = 0.1, value = 1),
	sliderInput("slider_alpha",  label = "Omega (ω)", min = 0, max = 1, step = 0.1, value = 0.1),
	sliderInput("slider_sigma",  label = "Sigma", min = 0, max = 1, step = 0.1, value = 0.1)
)

plot.inputs.bi <- box(width = 12,
									 sliderInput("slider_r1", label = "R1", min = 0, max = 1, step = 0.1, value = 0.1)
)

plot.params.bi <- box(width = 12,
									 sliderInput("slider_kappa",  label = "Kappa", min = 0, max = 1, step = 0.1, value = 1),
									 sliderInput("slider_gamma1", label = "Gamma", min = 0, max = 1, step = 0.1, value = 1),
									 sliderInput("slider_alpha",  label = "Omega (ω)", min = 0, max = 1, step = 0.1, value = 0.1),
									 sliderInput("slider_sigma",  label = "Sigma", min = 0, max = 1, step = 0.1, value = 0.1)
)


sidebar <- dashboardSidebar(
	sidebarMenu(
		id = "sidebarmenu",
		menuItem("About this App", tabName = "home",  icon = icon("home")),
		menuItem("Replication Model", tabName = "replication",  icon = icon("clone"), selected = T),
		menuItem("Bipartite Model", tabName = "bipartite", icon = icon("retweet")),
		menuItem("Tripartite Model", tabName = "tripartite", icon = icon("spinner")),

		# Replication Model
		conditionalPanel(
			condition = "input.sidebarmenu == 'replication'",
			plot.checkbox,
			actionButton("run_sims", label = "Run Simulations")
		),

		# Bipartite Model
		conditionalPanel(
			condition = "input.sidebarmenu == 'bipartite'",
			plot.checkbox,
			actionButton("run_sims2", label = "Run Simulations")
		),

		# Tripartite Model
		conditionalPanel(
			condition = "input.sidebarmenu == 'tripartite'",
			plot.checkbox
		)

	)
)

body <- dashboardBody(
	useShinyjs(),

# 	tags$head(
# 		tags$style(type="text/css",
# 							 "label.control-label, .selectize-control.single {
# 							 display: table-cell;
# 							 text-align: center;
# 							 vertical-align: middle;
# 							 }
# 							 label.control-label {
# 							 padding-right: 10px;
# 							 }
# 							 .form-group {
# 							 display: table-row;
# 							 }
# 							 .selectize-control.single div.item {
# 							 padding-right: 15px;
# 							 }")
#   ),

	tags$head(tags$style(
		HTML('
			.content-wrapper,
			.right-side {
				background-color: #ffffff;
			}
				 '))),

	tabItems(
		tabItem(tabName = "home",
						h2("About this project"),
						includeMarkdown("body.md"),
						actionButton(
							inputId='ab1', label="Fork us",
							icon = icon("github"),
							onclick ="window.open(
							'https://github.com/aldomann/open-data-project', '_blank')"
							)
		),

		# Replication Model
		tabItem(tabName = "replication",
						h2("Replication Model"),
						verbatimTextOutput("debugger"),
						fluidRow(
							column(width = 4,
								plot.inputs,
								plot.params
							),
							column(width = 8,
										 uiOutput('ui_plots')
							)
						)
		),

		# Bipartite Model
		tabItem(tabName = "bipartite",
						h2("Bipartite Model")
						# verbatimTextOutput("debugger"),
						# fluidRow(
						# 	column(width = 4,
						# 				 plot.inputs.bi,
						# 				 plot.params.bi
						# 	),
						# 	column(width = 8,
						# 				 uiOutput('ui_plots')
						# 	)
						# )

		),

		# Tripartite Model
		tabItem(tabName = "tripartite",
						h2("Tripartite Model")
						# dygraphOutput("compare.plot")
		)
	)
)

ui <- dashboardPage(
	skin = "green",
	header,
	sidebar,
	body
)
