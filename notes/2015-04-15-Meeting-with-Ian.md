# Meeting with Ian
_April 15, 2015_

Updates

* Turned conference work into a slide deck, practiced Quals 2 at group meeting. I don't have a Quals 2 committee yet, will follow up soon.

* Tips on the Quals 2 document?
	* Where does your work fit into the context of other work? Are you so focused that you only know your own work? If I were on the committee, I would see how well the student knows the optimization frameworks they are using.

	* Build the introduction to set the context. Liberally add references to relevant material.

	* Johanna will probably be on the committee. Others won't be particularly interested in instanton idea. They are more interested in:
		* how well do you understand the ideas? 
		* How well can you extend those ideas? 
		* How well can you present your body of knowledge? 
		* Can I go off on my own and make significant contributions?

	* Our paper lacked a substantial set of references. Given that we want to turn this into a journal version, there should be more references.

* Annual PhD progress report. I fill out my part, then send to Ian, right?
	* Yes. Fill out my part, send to Ian. Jon sent the PDF with the forms, Ian filled it out, printed it, signed and scanned it, then sent it back.

* PowerTech: Paper accepted? Hotel stay is covered 4 nights. One month until late registration fee kicks in.

* Line current calculation: Part 4 of IEEE 738 is most relevant.
	* Starting point is to say: if I have a line at steady state temperature and apply a step change to current, what does the temperature look like? The time constant is important. If it's 15 minutes, we don't need to look at time periods of 10 minutes or 2 hours. What is that time constant? That determines the sensible horizon to look at.
	* Once we know the horizon and time constants, we can look at what makes sense with Euler integration. If wind forecast comes at 5 minute basis, and time constants are 20 minutes, we don't care about 1 minute basis. But maybe we do care about five minutes, maybe 10 minutes is too slow.
	* We need to give greater weight to earlier time steps; uncertainty increases as we look further out.
	* Our temperature equation is not a linearization. Mads linearized and then used Euler integration. We should keep the nonlinear equation and use Euler integration on that.
	* Non-negativity constraints on wind generation. We really need these. If we get a solution with negative generation, we can truncate and re-run the problem. But actually we're only really concerned with cases where deviations are small. But what if forecasts themselves are small?
	* If Dan's model accommodates inequality constraints, we should have no problem.

* When are you getting back?
	* June 2. Other travel happens then, but Ian will be based in Ann Arbor then.

* Ian: the thermal model in the PowerTech paper is not correct. Work through that and send me what that model should look like? Get that to me soon. That's the fundamental part of the whole thing. Include details on time constants and line temperature. There's some offset that's not correct there.

* Quals 2: the PowerTech paper was light on results. Build that section of the report and presentation. The figure that shows the network -- I have no idea what that's showing.