//
// Notes.sc:
// ---------
//
// author:	Rainer Doemer
// last update:	01/16/02
//
// note:	this SpecC example demonstrates the syntax
//		of persistant annotation in SpecC code;
//		it has no other useful functionality


/*** general notes (attached to the design) ***/

note Author	= "Rainer Doemer";
note Update	= "Wed Jan 16 11:21:03 PST 2002";
note Date	= { 2002, 01, 16 };
note Comment	= "This is an example of how to attach notes to a design.\n"
		  "A note can be any constant expression and can be attached "
		  "to any symbol, label, or user-defined type in the code.";
note Comment2	= { "NEW:",
		    "Annotations now can be of composite type!",
		    "That is, the value of annotations can be structured "
		    "just like complex initializers for arrays or structs.",
		    "This change is effective as of SpecC LRM V2.0." };

note TicTacToe	= { { 1, 0, 1 },
		    { 0, 1, 0 },
		    { 1, 0, 1 } };


interface I
{
	void Write(int Value);

	/*** notes inside interface ***/
	note Name = "I";
	note Write.Complexity = 0;
	note Write.ArgType = "integer";
};


channel C implements I
{
	int	Storage;
	bool	Valid;

	void Write(int Value)
	{
	/*** notes inside member function ***/
	note Write.Name = "Write";
	note Value.Type = "int";

	l1: Storage = Value;
	l2: Valid = true;

	/*** notes at statements ***/
	note Storage.Color = "green";
	note l1.OpID = 42;
	note l2.OpID = 42+1;
	}

	/*** notes inside channel ***/
	note Name = "C";	// equivalent to: note C.Name = "C";
	note Storage.Bytes = sizeof(int);
};

behavior B(in int p1, I p2)
{
	int square(int a)
	{
	int r;
		{ int t;
		  /*** notes in compound statements ***/
		  note t.comment = "t is temporary only";
		  statement1: t = a * a;
		  note statement1.cost = 32 * 32;
		  statement2: r = t;
		 }
	statement3: return(r);
	}

	void main(void)
	{
	const int x = 5;
	int y;

	/*** notes inside function ***/
	note Name = "main";	// equivalent to: note main.Name = "main";
	note x.Magic = sizeof(int) * 16 / 2 + (-42)/6 + 7.0;

	label1: y = square(x);
	}

	enum { E_TEMPORARY, E_STANDARD, E_INLINED };

	/*** notes inside behavior ***/
	note B.Name = "B";
	note main.FullName = "B.main";
	note square.type = E_INLINED;
	note E_TEMPORARY.blablabla = "some note at enumerator constant";
};

/*** notes attached to interface/channel/behavior names ***/

note I.comment	= "This is an one-way interface.";
note C.cost	= 512.50;
note B.version	= 1.0;


/*** notes at aggregate types ***/

struct S
{
	int a, b;
	note a.bits	= sizeof(int) * 8;
	float f;
	note f.bits	= sizeof(float) * 8;
};

note S.comment	= "for demonstration only";
note S.bits	= sizeof(struct S) * 8;

typedef struct /* unnamed */
	{ float x, y, z; } Vector3D;
note Vector3D.bits = sizeof(Vector3D) * 8;


behavior Main	/* empty Main behavior */
{
	int main(void)
	{
	note Functionality = 0;
	return(0);
	}
};

// EOF
