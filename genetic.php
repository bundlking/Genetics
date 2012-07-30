<html>
<head>
<style>
	body {
		font-size:10px;
	}
	.arraybox {
		width:30%;
		float:left;
		border:solid 1px #3399FF;
		margin-left:10px;
	}
</style>    
</head>
<body>
<?php

class Algo {
	// define properties
	public $popsize;
	public $generation;
	public $maxiterations;
	public $eliterate;
	public $mutationrate;
	public $stop = false;
	public $target = array();
	public $genearray = array();
	public $parentarray = array();
	public $childarray = array();
	public $solution = array();
	
	// define functions
	public function initialise() {
		// initialise key variables for the algorithm
		$this->generation = 1;
		$this->popsize = 150;
		$this->maxiterations = 10000;
		$this->eliterate = 0.2;
		$this->mutationrate = 0.01;
		$this->target = str_split("The great Morpheus. We meet at last.");	
	}	

	// stopping condition test
	public function successtest() {
		foreach($this->genearray as $gene) {
			if($gene->fitness==count($this->target)) {
				$this->stop=true;
				$this->solution = $gene->val;
			}
		}
	}
	
	// stall condition test
}


class Gene {
	// define properties
	public $val = array();
	public $fitness;
	public $probability;
	public $probabilitysum;
	
	// define functions	
	// starter
	public function starter($algo) {
		
		while(count($this->val) < count($algo->target)) {
			// random character
			$arr = str_split('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz !@#$%^&.*()-+,_/?";:~`'); // get all the characters into an array
			shuffle($arr); // randomize the array
			$c = $arr[0]; // get the first (random) character out
			// add it to the value array of the gene
			$this->val[] = $c;
		}
		//print_r($this->val);
		//echo "<br />";		
	}

	public function uniform_crossover($parent1, $parent2) {
		for($x=0;$x<count($parent1->val);$x++) {
			// generate a random number between 0 and 1
			$random = rand(0,100)/100;			
			if($random<0.51) {
				$this->val[] = $parent1->val[$x];
			}
			else {
				$this->val[] = $parent2->val[$x];			
			}
		}
	}

	public function fitnesstest($algo) {
		$fitness = 0;
		
		for($x=0;$x<count($algo->target);$x++) {
			if($algo->target[$x]==$this->val[$x]) {
				// then a match on this letter - add to their fitness score
				$fitness++;				
			}
			else {
			}
		}		
		$this->fitness = $fitness;
	}

	public function mutate($algo) {
		// loop through every cell in the array, mutate if random number test passes
		for($x=0;$x<count($this->val);$x++) {
			// generate a random float between 0 and 100,000
			$random = rand(0,100000);
			// normalise to between 0 and 1
			$random = $random / 100000;
			// if $random < mutation % nominated for algorithm, then mutate
			if($random < $algo->mutationrate) {
				// select random new character to replace this value				
				$arr = str_split('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz !@#$%^&.*()-+,_/?";:~`'); // get all the characters into an array
				shuffle($arr); // randomize the array
				$c = $arr[0]; // get the first (random) character out
				// replace this value in the gene
				$this->val[$x] = $c;				
			}
		}
	}


	
}

// function used to sort by fitness
function fitSort( $a, $b ) {
    return $a->fitness == $b->fitness ? 0 : ( $a->fitness > $b->fitness ) ? 1 : -1;
}

$algo = new Algo;
$algo->initialise();

$count = 0;
while($count < $algo->popsize) {
	// create enough starting genes to get going
	$gene = new Gene;
	$gene->starter($algo);
	$gene->fitnesstest($algo);
	$algo->genearray[] = $gene;
	$count++;
}

while($algo->generation < $algo->maxiterations && $algo->stop==false) {
	// sort genearray by fitness (smallest -> largest)
	usort($algo->genearray, 'fitSort');
	
	$fitsum = 0;
		
	// ELITES
	// top solutions become children directly - calc from %, grab from the end of the genearray, add directly to childarray
	$numelites = $algo->eliterate * count($algo->genearray);
	$insertcount = count($algo->genearray) - 1;	
	while(count($algo->childarray)<$numelites) {
		$algo->childarray[] = $algo->genearray[$insertcount];
		$insertcount--;
	}
	
	// MUTANTS
	// mutate top solutions (elites) to offer speed to algorithm
	foreach($algo->childarray as $elite) {
		$mutant = new Gene;
		$mutant->val = $elite->val;
		$mutant->mutate($algo);
		// test fitness
		$mutant->fitnesstest($algo);
		$algo->childarray[] = $mutant;		
	}
	
	// ROULETTE SELECTION OF PARENTS
	// parents selection using weighted fitness - higher fitness more likely to be selected
	// http://www.cs.nott.ac.uk/~gxk/courses/g5baim/004ga/GA05-populationmod.html
	/* Roulette Wheel Selection : The idea behind the roulette wheel selection parent selection technique is that each individual is given a chance to become a parent in proportion to its fitness evaluation. It is called roulette wheel selection as the chances of selecting a parent can be seen as spinning a roulette wheel with the size of the slot for each parent being proportional to its fitness. Obviously those with the largest fitness (and slot sizes) have more chance of being chosen.
		Roulette wheel selection can be implemented as follows
		1. Sum the fitnesses of all the population members. Call this TF (total fitness).
		2. Generate a random number m, between 0 and TF.
		3. Return the first population member whose fitness added to the preceding population members is greater than or equal to m.
	*/
	foreach($algo->genearray as $g) {
		$fitsum = $fitsum + $g->fitness;
	}
	
	$probabilitysum = 0;
	foreach($algo->genearray as $g) {
		if($fitsum>0) {
			$g->probability = $g->fitness / $fitsum;
			$g->probabilitysum = $probabilitysum + $g->probability;
			$probabilitysum = $g->probabilitysum;
		}
	}
	
	// fill remainder of child array by selecting pairs of parents, using one point crossover for breeding and applying mutation to children
	while(count($algo->childarray)<count($algo->genearray)) {
		while(count($algo->parentarray)<2) {
			$selection = rand(0,100)/100;
			//echo "<h2>Selection: ".$selection."</h2>";
			$selected = false;
			for($x=0;$x<=count($algo->genearray)-1;$x++) {	
				//echo "test: ".$algo->genearray[$x]->probabilitysum."<br />";
				if($algo->genearray[$x]->probabilitysum > $selection && $selected==false) {
					//echo "SELECTED<br />";
					$algo->parentarray[]=$algo->genearray[$x];
					$selected=true;
				}
			}
		}
		
		// ONE POINT CROSSOVER TO BREED CHILDREN
		/*
		// 1. Select a random point for the crossover in the target
		$crossoverpoint = rand(0,count($algo->target)-1);
		
		$left = array();
		$right = array();
		
		// slice parent 'val' values
		foreach($algo->parentarray as $p) {
			//echo "PARENT: ".implode($p->val)."----------------------"; 
		
			$left[] = array_slice($p->val, 0, $crossoverpoint); 
			$right[] = array_slice($p->val, $crossoverpoint);
		}
		
		// Crossover to create children
		$child1 = new Gene;
		$child1->val = array_merge($left[0], $right[1]);
		// apply mutation
		$child1->mutate($algo);		
		// test fitness
		$child1->fitnesstest($algo);
		$algo->childarray[] = $child1;
	
		$child2 = new Gene;
		$child2->val = array_merge($left[1], $right[0]);
		// apply mutation
		$child2->mutate($algo);	
		// test fitness
		$child2->fitnesstest($algo);
		$algo->childarray[] = $child2;
		*/
		
		//echo "<br />Child 1: ".implode($child1->val)."----------------------";
		//echo "Child 2: ".implode($child2->val)."<br />";
		
		// UNIFORM CROSSOVER TO BREED CHILDREN
		// Every item in each child could come from either parent
		$child1 = new Gene;
		// apply crossover
		$child1->uniform_crossover($algo->parentarray[0], $algo->parentarray[1]);
		// apply mutation
		$child1->mutate($algo);		
		// test fitness
		$child1->fitnesstest($algo);
		$algo->childarray[] = $child1;
		
		$child2 = new Gene;
		// apply crossover
		$child2->uniform_crossover($algo->parentarray[0], $algo->parentarray[1]);
		// apply mutation
		$child2->mutate($algo);	
		// test fitness
		$child2->fitnesstest($algo);
		$algo->childarray[] = $child2;			
		
		// flush parent array for next iteration
		$algo->parentarray = array();
						
	}
	
	echo "<h1 style='clear:both;'>GENERATION: ".$algo->generation."</h1>";
	echo "<div class='arraybox' id='genearray'>";
	echo "<h2>Gene Array</h2>";
	//var_dump($algo->genearray);
	foreach($algo->genearray as $g) {
		echo implode($g->val)."<br />";
	}
	echo "</div>";	
	
	/*
	echo "<div class='arraybox' id='parentarray'>";
	echo "<h2>Parent Array</h2>";
	var_dump($algo->parentarray);
	echo "</div>";
	*/
	
	echo "<div class='arraybox' id='childarray'>";
	echo "<h2>Child Array</h2>";
	//var_dump($algo->childarray);
	foreach($algo->childarray as $g) {
		echo implode($g->val)."<br />";
	}	
	echo "</div>";
	
	$algo->generation++;
	
	// reset genearray
	$algo->genearray = array();
	
	// set genearray to the parent array
	$algo->genearray = $algo->childarray;
	
	// reset child array
	$algo->childarray = array();
	
	// test for success
	$algo->successtest();
	
}

if(count($algo->solution)>0) {
	echo "<h1>SOLUTION: ".implode($algo->solution)."</h1>";
}

?>

</body>
</html>